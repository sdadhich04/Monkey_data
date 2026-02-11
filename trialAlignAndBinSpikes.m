function [trial_spike_rate, time_bins] = ...
    trialAlignAndBinSpikes(spike_file, align_event_code, time_before, time_after, bin_width)

% [trial_spike_rate, aligned_spike_times, aligned_spike_labels, time_bins] = ...
%    trialAlignAndBinSpikes(spike_file, align_event_code, time_before, time_after, bin_width)
%
%Function to align spikes to trial-event-times and compute trial-aligned
%firing rates. 
%A. Orsborn (created 12/27/19, updated 1/2/20)
%
%inputs: spike_file - string with full path to spike file to load
%                     (containing 'sig' variables, 'EVENTS', and 'EVENT_TIMES')
%        align_event_code - scalar, value of event-code to align spikes to.
%        time_before      - time prior to alignment event to load for a
%                           given trial (in seconds)
%        time_after      - time after alignment event to load for a
%                           given trial (in seconds)
%        bin_width       - size of bin to use for estimating spike rate (in seconds)
%
%outputs: trial_spike_rate    - matrix [#trials x time x #units] of estimated spike-rates for each trial/unit
%         time_bins           - vector (#time points x 1) of the center of
%                               the time-bin over which firing rate is estimated. (in seconds)
%


%load the spike_times, EVENTS and EVENT_TIMES from the spike-file
load(spike_file, 'spike_times', 'EVENTS', 'EVENT_TIMES')

%get the list of all neurons ('sig*' variables)
num_units = length(spike_times);


%%%%%% trial-alignment for each spike

%find the alignment-event times (i.e. time-stamp when EVENTS = align_event_code)
align_time = EVENT_TIMES(EVENTS == align_event_code);
num_trials = length(align_time);

%trial-align the spike times for each unit
%hint: trialAlignSpikes.m from lab2

% define bin edges and centers
edges = (-time_before):bin_width:(time_after);
time_bins = edges(1:end-1) + bin_width/2;   % bin centers
num_bins = length(time_bins);

% preallocate output: trials x time x units
trial_spike_rate = zeros(num_trials, num_bins, num_units);

% loop through units
for iU = 1:num_units

    st = spike_times{iU};  % spike times for this unit (in seconds)

    % loop through trials
    for iT = 1:num_trials

        % align spikes to this trial event
        aligned = st - align_time(iT);

        % keep only spikes in the window
        aligned = aligned(aligned >= -time_before & aligned <= time_after);

        % bin spikes (counts per bin)
        counts = histcounts(aligned, edges);

        % convert counts to rate (spikes/sec)
        trial_spike_rate(iT,:,iU) = counts ./ bin_width;
    end
end



%now bin your trial-aligned spikes
%you've already written code to do this. 

