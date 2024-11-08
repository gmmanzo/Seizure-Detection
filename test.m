clear 
clc 
% Load EDF file using biosig
file = uigetfile('*.edf.seizures'); % Opens a dialog box in current folder to select file for use
[data, info] = sload(file);

% Display basic information about the EDF file
disp('File Information:');
disp(info);

time = (0: length(data)-1)/info.SampleRate;
numSignals = size(data, 2);


%creates one big figure of each channel 
figure;
for i = 1:numSignals
    subplot(numSignals, 1, i);
    plot(time, data(:,i));
    title(['Signal ', num2str(i), ': ', info.Label{i}]);
    xlabel('Time (s)');
    ylabel('Amplitude');
end
%% 


% Plot EEG signals on single graph 
time = (0: length(data)-1)/info.SampleRate; % Time vector based on sampling rate
figure; % Opens a figure window
plot(time, data);
xlabel('Time (s)');
ylabel('Amplitude (uV)');
title('EEG Signals');
legend(info.Label); % Adds labels as legend for each signal

% Bandpass filter for Gamma waves
fs = info.SampleRate; 
[b, a] = butter(4, [30 50]/(fs/2), 'bandpass') % Fourth order butterworth filter for a steeper cutoff, normalized by nyquist frequency
filter_EEG = filtfilt(b, a, data); % Zero phase filtering

figure; 
plot(time, filter_EEG);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('EEG data with Gamma bandpass filter');
legend(info.Label);

% Wavelet transform filtering
filter_EEG_wave = wdenoise(data,5); % Level 5 wavelet
figure; 
plot(time, filter_EEG_wave);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('EEG data with wavelet filter');
legend(info.Label);

%% 


% Plot the signals
% numSignals = size(record, 2);  % Number of signals in the EDF file
% time = (0:size(record, 1)-1) / 256;  % Time vector
% %time = (0:size(record, 1)-1) / header.SampleRate(1);  % Time vector

%Signal filtering data
% Design a low-pass filter (e.g., 50 Hz cutoff)
Fc = 1;                  % Cutoff frequency
Fs = 256;
[b, a] = butter(4, Fc/(Fs/2), 'low');  % 4th order low-pass Butterworth filter

% Apply the filter to the signal
filtered_signal = filtfilt(b, a, record(:,1));

% Plot the original and filtered signal
figure;
plot(time, record(:,1)); hold on;
plot(time, filtered_signal, 'r');  % Filtered signal in red
xlabel('Time (s)');
ylabel('Amplitude');
title('Original vs Filtered Signal');
legend('Original', 'Filtered');

% Create separate figures for each signal
    for i = 1:numSignals
        figure;  % Create a new figure for each signal
        plot(time, record(:,i));
        
        if isfield(header, 'Label') && iscell(header.Label) && numel(header.Label) >= i
            title(['Signal ', num2str(i), ': ', header.Label{i}]);
        else
            title(['Signal ', num2str(i)]);
        end
        
        xlabel('Time (s)');
        ylabel('Amplitude');
    end


% Adjust the plot layout
sgtitle('EDF File Signals');


%creates one big figure of each channel 
figure;
for i = 1:numSignals
    subplot(numSignals, 1, i);
    plot(time, record(:,i));
    title(['Signal ', num2str(i), ': ', header.Label{i}]);
    xlabel('Time (s)');
    ylabel('Amplitude');
end
