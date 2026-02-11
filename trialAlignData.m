function [trial_data] = trialAlignData(data, trial_event_times, time_before, trial_segment_time, Fs)


%trialAlignData   (A. Orsborn, created: 4-22-10, updated 12/26/19)
%[trial_data] = trialAlignData(data, trial_event_times, time_before, trial_segment_time, Fs)
%aligns a data time-series (data) to the trial_event_times. Returns
%time-series data segments that are trial_segment_time long for each trial that begin time_before prior to
%the trial-event
%
%input:
%       data               - data to trial sort. (time x n) matrix, where n is # variables. 
%       trial_event_times  - time of trial events to align to (#trials x 1). Time units must be
%                            consistent with sampling rate (e.g. if event-times in seconds, Fs must be in Hz)
%       time_before        - time before event-code occurence to include in data (in units consistent with Fs)
%       trial_segment_time - total length of trial-aligned data segments (in units consistent with Fs)
%       Fs    - sampling rate
%output:
%       trial_data - Trial-aligneddata (trials x time x n)


nTr  = length(trial_event_times); %# of trials


%initialize trial_data
trial_data = zeros(nTr, trial_segment_time*Fs, size(data,2));


%loop through trials
for iT=1:nTr
    
    %select time-indices for this trial
    L1 = floor(trial_event_times(iT)* Fs) - floor(time_before*Fs);
    L2 = L1 + floor((trial_segment_time)*Fs)- 1;

    %confirm within bounds of data (only happens for incorrectly formatted data)
    if L1 > 0 & L2 > 0 & L2 <= size(data,1) %#ok<AND2>
        
        if L2 <= size(data,1)
            trial_data(iT,:,:) = data(L1:L2,:);
        end
        
    else
        disp(strcat('Error in time indexing. Skipping trial ', num2str(iT)))
    end
end
