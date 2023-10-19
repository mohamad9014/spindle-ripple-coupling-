function [spindles] = find_spindles(data, fs, freq_range, threshold, duration_range, merge_time, ch_name)
% FIND_SPINDLES  Find spindles.
%   SPINDLES = FIND_SPINDLES(DATA,FS,FREQ_RANGE,THRESHOLD,DURATION_RANGE,MERGE_TIME) finds spindles in DATA that sampled with frequency FS [Hz].
%   FREQ_RANGE [Hz] is the frequency range of spindles.
%   THRESHOLD is a factor of the standard deviation of filtered DATA.
%   Events merge when they are too close (less than MERGE_TIME [s]).
%   Events discard when they are too long or short (DURATION_RANGE [s]).
%
%   SPINDLES = FIND_SPINDLES(DATA,FS,FREQ_RANGE,THRESHOLD,DURATION_RANGE,MERGE_TIME,CH_NAME) If CH_NAME is set, the algorithm will be shown in the plot.
%
%   SPINDLES is a struct that contains FILTERED_DATA, SPINDLE_DATA, START,
%   STOP, and PEAKS_IND.
%
%   EXAMPLE:
%       spindles = find_spindles(data, 1000, [7 15], 2.5, [0.5 3], 0.5, 'Thalamus');
%
%   Reference:
%       Alizadeh, Z., Azimi, A., Ghorbani, M. (2021). Enhancement of 
%       hippocampal-thalamocortical temporal coordination during 
%       slow-frequency long-duration anterior thalamic spindles. 
%       BioRxiv, 2021.10.03.462943. 
%       https://doi.org/10.1101/2021.10.03.462943



time = (1:length(data))/fs;

%% Filter
% The data is first bandpass filtering in the spindle frequency range
% using FIR filters from the EEGLAB toolbox. (filter order corresponds to 3
% cycles of the low frequency cut off)
f_data = eegfilt(data, fs, freq_range(1), freq_range(2), 0, fix(fs/freq_range(1))*3, 0, 'fir1');

%% Envelope
% Hilbert transform is used to compute instantaneous amplitude.
e_data = abs(hilbert(f_data));

%% Smooth
% Enveloped data is smoothing using a 300 ms Gaussian window.
s_data = smoothdata(e_data, 'gaussian', round(0.3*fs));

%% Threshold
sd_data = std(f_data);
th = threshold*sd_data;

%% Find start/stop as indices of threshold crossings
crossings = diff(s_data > th);

start = find(crossings>0);
stop = find(crossings<0);

if time(stop(1)) < time(start(1))
    stop(1) = [];
end
if time(start(end)) > time(stop(end))
    start(end) = [];
end

% Merge events when they are too close
while true
	tooClose = find(time(start(2:end)) - time(stop(1:end-1)) < merge_time);
	if isempty(tooClose), break; end
	start(tooClose+1) = [];
	stop(tooClose) = [];
end

% Discard events that are too long or too short
duration = time(stop) - time(start);
bad = duration < duration_range(1) | duration > duration_range(2);
start(bad) = [];
stop(bad) = [];

spindle = zeros(size(data));
peaks_ind = zeros(size(start));
for i=1:length(start)
    spindle(start(i):stop(i)) = 1;
    
    [~,ind] = max(f_data(start(i):stop(i)));
    peaks_ind(i) = start(i)+ind-1;
end

spindles.filtered_data = f_data; 
spindles.spindle_data = spindle;
spindles.start = start;
spindles.stop = stop;
spindles.peaks_ind = peaks_ind;

%% plot
if nargin > 6
    figure
    plot(time, f_data)
    hold on
    plot(time, s_data)
    plot(time, th*ones(length(data),1))
    plot(time, spindle*max(s_data))
    legend(ch_name, 'Smoothed envelope', 'Threshold', 'Spindle')
    xlabel('Time [s]')
end