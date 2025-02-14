import tensor_flow as tf
import numpy as np

def extract_features(windows):
    features = {}

    # Time-domain features
    features['mean'] = np.mean(windows)
    features['std'] = np.std(windows)
    features['skew'] = skew(windows)
    features['kurtosis'] = kurtosis(windows)

    # Frequency-domain features (Power Spectral Density)
    f, Pxx = welch(windows, fs=512, nperseg=512)
    features['power_band'] = np.sum(Pxx)  # Total power

    return features
     

def edf_file_extractor(file_path):
  """
  Opens an EDF file and returns its contents in the form of a Data Frame
  """
  # Load all data
  import pyedflib # Used for edf files
  import pandas as pd

  # Load the EDF file
  edf = pyedflib.EdfReader(file_path, check_file_size=False)

  # Extract signal labels
  signal_labels = edf.getSignalLabels()

  # Extract signal data
  signals = []
  for i in range(edf.signals_in_file):
      signals.append(edf.readSignal(i))

  # Create a DataFrame
  df = pd.DataFrame(signals).transpose()
  df.columns = signal_labels

  # Close the EDF file
  edf.close()

  return df
     

def seizure_time_parser(txt_file_path):
    """
    Takes in the txt file in the dataset and determines
    the times of the seizures (if they have any) for all
    files listed in the txt file.

    Returns a dictionary in the following format:
      {Filename: [[Seizure start time, seizure end time], ...]}
    """
    if not txt_file_path.endswith('.txt'):
        print("File is not a txt file")
        return {}

    seizure_dict = {}

    # Read non-empty stripped lines
    with open(txt_file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 0
    while i < len(lines):
        line = lines[i]
        # Look for a file entry
        if line.startswith("File Name:"):
            # Extract the file name
            filename = line.split(":", 1)[1].strip()
            # Advance until we find the "Number of Seizures in File:" line.
            i += 1
            while i < len(lines) and not lines[i].startswith("Number of Seizures in File:"):
                i += 1

            if i < len(lines):
                count_line = lines[i]
                count_str = count_line.split(":", 1)[1].strip()
                try:
                    seizure_count = int(count_str)
                except ValueError:
                    seizure_count = 0

                # Only process if at least one seizure exists.
                if seizure_count > 0:
                    seizure_intervals = []
                    # Process each seizure.
                    for j in range(seizure_count):
                        # The next line should indicate the seizure start time.
                        i += 1
                        if i < len(lines) and "Start Time:" in lines[i]:
                            start_line = lines[i]
                            # Split on ":" and take the first token of the remainder to get the numeric part.
                            parts = start_line.split(":", 1)
                            if len(parts) > 1:
                                start_time_str = parts[1].strip().split()[0]
                                try:
                                    seizure_start = int(start_time_str)
                                except ValueError:
                                    seizure_start = None
                            else:
                                seizure_start = None
                        else:
                            seizure_start = None

                        # The following line should indicate the seizure end time.
                        i += 1
                        if i < len(lines) and "End Time:" in lines[i]:
                            end_line = lines[i]
                            parts = end_line.split(":", 1)
                            if len(parts) > 1:
                                end_time_str = parts[1].strip().split()[0]
                                try:
                                    seizure_end = int(end_time_str)
                                except ValueError:
                                    seizure_end = None
                            else:
                                seizure_end = None
                        else:
                            seizure_end = None

                        # Only add the seizure interval if both times were successfully parsed.
                        if seizure_start is not None and seizure_end is not None:
                            seizure_intervals.append([seizure_start, seizure_end])
                    # Store the seizure intervals for the file.
                    seizure_dict[filename] = seizure_intervals
        i += 1

    return seizure_dict
     

# Power Spectral Density Version of Data:
def psd(node, plot=False, name="____",find_max=True):
  """
  Takes in a node, which refers to the specific EEG signal
  Calculates the Power Spectral Density of a given EEG signal.
  Returns raw data as f, Power.
  F represents Frequency
  Power is Power
  """
  from scipy.signal import welch
  import matplotlib.pyplot as plt

  f, Power = welch(node, fs=512, nperseg=512)

  if plot & find_max:
    plt.figure(figsize = (12,7))
    plt.plot(f, Power)
    plt.title(f"Power Spectral Density of {name}")
    plt.xlabel("Frequency(Hz)")
    plt.xlim(0,20)
    plt.plot(f[pd.DataFrame.idxmax(pd.DataFrame(Power))],max(Power),'ro')
  elif plot:
    plt.figure(figsize = (12,7))
    plt.plot(f, Power)
    plt.title(f"Power Spectral Density of {name}")
    plt.xlabel("Frequency(Hz)")
  return f, Power
     

# 1D Discrete Fast Fourier Transform Version of data:
def fft(node, plot=False, name="_____"):
  """
  Takes in a node, which refers to the specific EEG signal
  Calculates the Fast Fourier Transform of a given EEG signal.
  Returns raw fft data
  """
  from scipy.fft import fft, fftfreq, rfft
  import numpy as np
  import matplotlib.pyplot as plt
  import tensorflow as tf

  fft_data = tf.signal.rfft(tf.cast(EEG_Fp1, tf.float32)) # Need to convert the data type into a 64 bit to run the FFT function
  if plot:
    plt.figure(figsize = (12,7))
    plt.plot(fft_data)
    plt.xlim(0,20)
    plt.grid()
    plt.xlabel('Frequency')
    plt.ylabel('Amplitude')
    plt.title(f"1D FFT of {name} data")

  return fft_data
     

def preprocessing_EEG_signal(node, plot=False, name="______"):
  """
  Take in a node, which refers to the speicific EEG signal.
  Filters the data.
  """
