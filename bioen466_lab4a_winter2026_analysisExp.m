%%%%% This script provides a guided outline for orienting yourself to the
%%%%% dataset we will be analyzing for labs 4 and 5

%%%%% Written by A.L. Orsborn, v191223
%%%%%
%%%%%
%%%%% All lines where you have to fill in information is tagged with a comment including "FILLIN". Use this flag to find everything you need to modify.
%%%%% all figures that need to be included in comprehension questions are
%%%%% flagged with %INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS

%% data metadata 

% Define some basic things to make it easy to find data files.
% We will want to take advantage of systematic naming structure in our data files.

dataDir = 'C:\Users\metal\Downloads\lab_ee_466\lab4a_ee_466\data'; %FILLIN with the path to where your data is stored

codeDir = 'C:\Users\metal\Downloads\lab_ee_466\lab4a_ee_466'; %FILLIN with the path to where your provided functions for this lab are stored

%we now add all the functions in our code directory to the matlab path
cd(codeDir)
addpath(genpath(codeDir))


%now define which specific file you will analyze
fileBase = 'jeev050712a_bioe466_winter2020';

spikes_file_tag = '_spikes.mat';
kinematics_file_tag = '_kinematics.mat';
lfp_file_tag    = '_lfp.mat';


%define some constants that help with data interpretation
GOCUE_CODE      = 5; %event-code flagging when the go-cue comes on
LEAVE_CENTER_CODE= 6; %event-code flagging when hand leaves center
REWARD_CODE     = 9; %event-code flagging when reward turns on
AT_TARGET_CODE    = 7; %event-code flagging when the animal is at the target
TIMEOUT_CODE      = 8; %event-code flagging when a reach time-out (didn't get to target fast enough) error occurred

TARGET_OFFSET   = 63; %offset to subtract from event-codes for target direction 

definePlotColors; %script to set-up figure aesthetics. 

%% We will now load kinematic data to visualize and manipulate it. 

%%%%%%%%%%% A: load the kinematic data file and visualize basics of the data. 

load(fullfile(dataDir, [fileBase kinematics_file_tag])) %FILLIN with data file pointer (path + name)

%explore the data you loaded
%reminder, the hand kinematics are in cartesian coordinates 
%(X,y) position (cm), velocity (cm/s), & acceleration (cm/s^2) in each dimension (time x 6) 
whos hand_kinematics hand_kinematics_subsampled FS_kinematics FS_kinematics_sub EVENTS EVENT_TIMES


%define a time-vector for the hand kinematics (we'll just work with the sub-sampled kinematics for now)
time_kinematics_sub = (0:size(hand_kinematics_subsampled,1)-1)./FS_kinematics_sub; %FILLIN

%plot the hand position in cartesian space for time range t = [0 180],   (i.e. the first 3 minutes);
figure
plotinds = time_kinematics_sub < 180;
plot(hand_kinematics_subsampled(plotinds,1), hand_kinematics_subsampled(plotinds,2),"r") %FILLIN
xlabel('X position (cm)') %FILLIN, include units
ylabel('Y position (cm)') %FILLIN, include units
title('Hand position (first 3 minutes)')  %FILLIN
axis equal

%plot the hand speed over time
%compute speed
speed = sqrt(hand_kinematics_subsampled(:,3).^2 + hand_kinematics_subsampled(:,4).^2); %FILLIN
figure
plot(time_kinematics_sub, speed,"r")
set(gca, 'xlim', [5 15]) %zoom in some to see structure
xlabel('Time (s)') %FILLIN, include units
ylabel('Speed (cm/s)') %FILLIN, include units
title('Hand speed over time')  %FILLIN

%%%%%%%%%%% B: Now analyze data with respect to task-structure:
%
%we will now "sort" the data based on the trial-structure within the data.
%To do this, we need to use the 'EVENTS' and 'EVENT_TIMES'. 
%EVENTS = a list of integer 'codes', which convey meaning about what happened in the task. 
%EVENT_TIMES = a list of time-stamps for when each 'event' in EVENTS happened in the file. (in units of seconds)
%
%The lab manual defines all of the relevant event-codes and their meanings.
%We provide a function 'trialAlignData' that performs the trial-alignment
%operations. Look at the 'help' for this function and use it for the
%following steps. 
%

%Let's try aligning the hand speed to when the hand leaves the center target:
%1. find all the EVENT_TIMES when EVENTS is the code of interest. 
align_time = EVENT_TIMES(EVENTS == LEAVE_CENTER_CODE); %FILLIN

%2. use trialAlignData to align the data. Load segments that are 3 seconds
%long, and include 1.5s of data before the event.
time_before = 1.5; %FILLIN
segment_length = 3.0; %FILLIN
trial_speed_lc = trialAlignData(speed, align_time, time_before, segment_length, FS_kinematics_sub); %FILLIN

%3. plot your trial-aligned speed against time (mean + traces of single
%trials). Plot a sub-set of trials (nPlot) to make things easier to
%visualize
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
trial_time_kin = ((0:(segment_length*FS_kinematics_sub)-1)./FS_kinematics_sub) - time_before; %FILLIN
figure
nPlot = 75;
plot(trial_time_kin, trial_speed_lc(1:nPlot,:)'); %individual trials
hold on
plot(trial_time_kin, mean(trial_speed_lc(1:nPlot,:),1), 'k', 'linewidth', 3) %FILLIN. plot mean -- fill in data
plot([0 0], [0 max(speed)], 'k--'); %plots a black line at 0 to visualize alignment
xlabel('Time from leave-center (s)') %FILLIN, include units
ylabel('Speed (cm/s)') %FILLIN, include units
title('Hand speed aligned to leave center') %include alignment

%now also generate another version of this figure, but aligning speed to
%the go-cue:
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
%hint - you can repeat the code above but change which code you align the
%data to. Just watch your variable names to keep track of which alignment
%you're using (i.e. don't overwrite trial_speed_lc, make a new variable.)
align_time_gc = EVENT_TIMES(EVENTS==GOCUE_CODE);
trial_speed_gc = trialAlignData(speed, align_time_gc, time_before, segment_length, FS_kinematics_sub);
figure
nPlot = 75;
plot(trial_time_kin, trial_speed_gc(1:nPlot,:)'); hold on
plot(trial_time_kin, squeeze(mean(trial_speed_gc,1)), 'k', 'linewidth', 3);
plot([0 0], [0 max(speed)], 'k--');
xlabel('Time from GO CUE (s)'); ylabel('Hand speed (cm/s)');
title('Hand speed aligned to GO CUE');


%%%%%%%%%%% C: Further divide data by sorting trials into types:

%As we saw in the hand trajectories plot, the data also has additional
%structure because the reach target position varies across each reach. 
%we're now going to use the event-data to split up the trials into groups
%for each reach direction. 
%To do this, we will "trial-align" the event codes as well to see and use
%the trial's structure 
%We provide a function 'trialAlignEvents' that performs the trial-alignment
%operations for event-data. Look at the 'help' for this function and use it for the
%following steps. It's the same idea as the time-series data alignment, but
%for event-codes. 
%

%trial-sort the events to leaving the center, loading 3 events before and
%after
[trial_events, trial_event_times] = trialAlignEvents(EVENTS, EVENT_TIMES, LEAVE_CENTER_CODE, 3, 3); %FILLIN

%look at the first 10 trials:
trial_events(1:10,:)

%You should see a clear structure in the events:
% [# between 64- 71]  15 GO_CUE_CODE LEAVE_CENTER_CODE ENTER_TARGET_CODE 51  REWARD_CODE
% the first column (# between 64 - 71) tells us which target (1 - 8) is
% being presented. They are offset by TARGET_OFFSET = 63. 
%
%based on this information, create a variable that tells you the reach
%direction of each trial.
trial_reach_direction = trial_events(:,1) - TARGET_OFFSET; %FILLIN


%now plot the mean cartesian hand trajectory (across trials), separated by
%the reach direction:
%1. trial-align the hand kinematics
[trial_kinematics_subsampled] = trialAlignData(hand_kinematics_subsampled(:,1:2), align_time, time_before, segment_length, FS_kinematics_sub); %FILLIN

%2. loop through all of the directions and compute the average across trials
%for each direction and plot
direction_list = unique(trial_reach_direction); %list of all targets
num_directions = length(direction_list);        %# of targets

mean_position_by_dir = nan(num_directions, length(trial_time_kin), 2); %initialize data matrix (#dirs x time x 2)
figure
hold on
for iD=1:num_directions
    
    mean_position_by_dir(iD,:,:) = squeeze(mean(trial_kinematics_subsampled(trial_reach_direction==direction_list(iD),:,:),1)); %FILLIN
    plot(squeeze(mean_position_by_dir(iD,:,1)), squeeze(mean_position_by_dir(iD,:,2)), 'linewidth', 2, 'color', plotColors(direction_list(iD),:)); %FILLIN
end
xlabel('X position (cm)') %FILLIN, include units
ylabel('Y position (cm)') %FILLIN, include units
title('Mean hand trajectory by reach direction')  %FILLIN

%Now adjust the plotting settings in this code to eliminate the extra
%wiggles we see at the center target. 
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure; hold on
post0 = trial_time_kin >= 0; %only plot after leaving center to avoid center-target wiggles
for iD=1:num_directions
    plot(squeeze(mean_position_by_dir(iD,post0,1)), squeeze(mean_position_by_dir(iD,post0,2)), 'linewidth', 2, 'color', plotColors(direction_list(iD),:));
end
xlabel('X position (cm)') %FILLIN, include units
ylabel('Y position (cm)') %FILLIN, include units
title('Mean hand trajectory by reach direction(After centre target)')  %FILLIN

%% we will now load and visualize neural data (spiking activity of neurons)

%%%%%%%%%%% D: Load and get familiar with neural spiking data

%load the spiking data file
load(fullfile(dataDir, [fileBase, spikes_file_tag])) %FILLIN with data file pointer (path + name)

%explore the data you loaded
%reminder, the spiking data is stored as detected spike-times for each multi-unit/neuron.
%Each multi-unit/neuron is saved in a cell within spike_times {#units x 1}.
%Each entry in spike_times is a vector (# spikes x 1) with time-stamps (in seconds)
%of each detected spike 
whos spike_times EVENTS EVENT_TIMES


%%%%%%%%%%% E: analyze spiking data with respect to task structure

%we will now align spikes to trial events (similar to the kinematics)
%and then turn the spike-times into spike rates ("binning"). 
%This analysis is the same as what we did in lab #2. 

%We have provided a function trialAlignAndBinSpikes, which is partially
%completed. Fill this in and then use it to trial-align and bin your spikes. 
%load 1.5 seconds before and after the alignment-event. Start with a
%bin-width of 100ms (we will vary this)
spike_file = fullfile(dataDir, [fileBase spikes_file_tag]); %FILLIN
time_before = 1.5; %FILLIN, in seconds
time_after = 1.5;  %FILLIN, in seconds
bin_width  = 0.1;   %FILLIN, in seconds
[trial_spike_rate, trial_time_spikes] = ...
    trialAlignAndBinSpikes(spike_file, LEAVE_CENTER_CODE, time_before, time_after, bin_width);


%now try plotting the trial-averaged PSTH (peri-stimulus time-histogram)
%for a few units. We will use some example units in the data for
%illustration purposes. We'll use all the units later. 
plot_units = [2 17 27 34];

%to help visualize how neural activity relates to kinematics, 
%let's also have a subplot of the hand speed 
figure
subplot(2, 1, 1)
hold on
plot(trial_time_spikes, squeeze(mean(trial_spike_rate(:,:,plot_units),1)), 'linewidth', 2) %FILLIN. plot trial-averaged PSTH for select neurons vs time
plot([0 0], [0 22], 'k--'); %plot a black line at 0 to visualize alignment
set(gca, 'ylim', [0 22])
xlabel('Time from leave-center (s)') %FILLIN, include units
ylabel('Firing rate (sp/s)') %FILLIN, include units
title('Trial-averaged PSTH (aligned to leave center)')  %FILLIN

subplot(2,1,2)
hold on
plot(trial_time_kin, squeeze(mean(trial_speed_lc,1)), 'k', 'linewidth', 3) %FILLIN. plot mean hand speed vs. time
plot([0 0], [0 max(mean(trial_speed_lc,1))], 'k--'); %plot a black line at 0 to visualize alignment
xlabel('Time from leave-center (s)') %FILLIN, include units
ylabel('Speed (cm/s)') %FILLIN, include units
title('Mean hand speed (aligned to leave center)')  %FILLIN


%%%%%%%%%%% F: analyze spiking data across trial types

%now let's try doing the same thing for the neurons as we did for the
%kinematics--splitting the trials up into which direction the animal is reaching. 
%Pick 2 example units (from the list plot_units above) and plot the average PSTH separated by reach
%direction. To make it easier to visualize, only plot the odd reaction
%directions (i.e. targets 1, 3, 5, and 7)

%INCLUDE THESE FIGURES IN COMPREHENSION QUESTIONS
plot_units = [2 34];


for iU=1:length(plot_units) %loop through units to plot
    
    figure
    hold on
    
    %loop through directions (odd only) to get mean and plot
    %FILLIN
    odd_dirs = direction_list(mod(direction_list,2)==1); 
    for iD=1:length(odd_dirs)
        d=odd_dirs(iD);
        tr_idx=(trial_reach_direction==d);
        mean_psth=squeeze(mean(trial_spike_rate(tr_idx,:,plot_units(iU)),1));
        plot(trial_time_spikes, mean_psth);
    end
    %end loop
    
    %label your figure
    xlabel('Time from leave-center (s)') %FILLIN, include units
    ylabel('Firing rate (sp/s)') %FILLIN, include units
    title(['Mean PSTH, split by reach target. Unit', num2str(plot_units(iU))])
end %loop units


%try re-generating these plots varying the bin_width to understand the role
%of binning in our estimates of neuron firing. 



%%%%%%%%%%% G: Compute "summary metrics" of spiking relationship to
%%%%%%%%%%% movement (direction-tuning)

%Finally, let's try making a plot to summarize how neuron firing rates
%relate to reach direction. We'll specifically plot the "direction
%tuning"--average firing rate in a time window vs. movement direction. 

%1. pick a time-window of interest. Let's start with the time right around
%movement onset: [-.1 .1]
tuning_window = [-0.1 0.1]; %FILLIN [startTime endTime] vector
t_index = (trial_time_spikes>=tuning_window(1)) & (trial_time_spikes<=tuning_window(2)); %FILLIN. Create a logical vector that is 1 for all times within the tuning window and 0 otherwise. 

%2. Compute the mean firing rate in this time-window 
%hint: you should have a [# trials x neurons] matrix
mean_rate_window = squeeze(mean(trial_spike_rate(:,t_index,:),2)); %FILLIN

%3. compute the mean & standard ERROR across reach directions
mean_rate_window_by_dir = nan(num_directions, size(mean_rate_window,2)); %FILLIN. initialize your matrices
ste_rate_window_by_dir = nan(num_directions, size(mean_rate_window,2));    %FILLIN.

%loop through directions to calculate
%FILLIN
for iD=1:num_directions
    dir = direction_list(iD);
    tr_idx = (trial_reach_direction==dir);
    mean_rate_window_by_dir(iD,:) = mean(mean_rate_window(tr_idx,:),1);
    ste_rate_window_by_dir(iD,:) = std(mean_rate_window(tr_idx,:),0,1)./sqrt(sum(tr_idx));
end %FILLIN


%4. plot the firing rate vs. reach direction for select units w/
%standard-error errorbars
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
plot_units = [2 17 27 34];

figure
errorbar(direction_list, mean_rate_window_by_dir(:,plot_units), ste_rate_window_by_dir(:,plot_units)) %FILLIN
xlabel('Reach direction (target #)') %FILLIN
ylabel('Mean firing rate in window (sp/s)') %FILLIN
title(['Direction tuning, select units. Time window = ', num2str(tuning_window(1)) '-', num2str(tuning_window(2)) 'S'])


% try re-running this analysis with varying time-windows. 


%% Finally, we will look at local field potential (LFP) data recorded from the same electrodes as the unit activity. 

%%%%%%%%%%% H: Load and get familiar with neural spiking data

%load the LFP data file
load(fullfile(dataDir, [fileBase, lfp_file_tag])); %FILLIN with data file pointer (path + name)

%explore the data you loaded
%reminder, the LFP data is saved as a matrix (# electrodes x time) with voltages (units: mV) 
whos LFP EVENTS EVENT_TIMES FS_lfp


%define a time-vector for the LFP
time_lfp = (0:size(LFP,2)-1)./FS_lfp;
%time_kinematics_sub = ; %FILLIN


%plot the LFP for a few electrodes for time range t = [0 5], (i.e. the first 10 seconds);
plot_electrodes = [1 5 10];
figure
plotinds = ismember(1:size(LFP,1), plot_electrodes); %FILLIN. logical vector that is 1 for the electrodes you want to plot and 0 for the ones you don't want to plot. 
tinds = time_lfp < 5;
plot(time_lfp(tinds), LFP(plotinds, tinds)'); %FILLIN. Plot the LFP data specified. 
xlabel('Time (s)') %FILLIN. include units
ylabel('LFP (mV)') %FILLIN. 
title('LFP, select electrodes (first 10 s)')  %FILLIN. 

%try zooming in a bit. What do you generally notice across the electrodes?


%%%%%%%%%%% I: analyze LFP in relationship to task structure

%Now let's trial-align the LFP to see if there are consistent relationships
%with behavior. 

%1. trial-align to leaving the center (load the same length of data as we
%have for kinematics and spikes)
%FILLIN
lfp_time_before = 1.5;
lfp_segment_length = 3;
trial_lfp = trialAlignData(LFP', align_time, lfp_time_before, lfp_segment_length, FS_lfp);
trial_time_lfp = ((0:(lfp_segment_length*FS_lfp)-1)./FS_lfp) - lfp_time_before; %FILLIN

%2. pick a few electrodes (suggestion: 1 and 5) and plot the LFP vs. time averaged across all
%trials. Reminder to have the x-axis be time in seconds. 
%FILLIN
plot_electrodes2 = [1 5];
figure;
hold on;
for e=plot_electrodes2
    plot(trial_time_lfp, squeeze(mean(trial_lfp(:,:,e),1)));
end

plot([0 0], [min(ylim) max(ylim)], 'k--');
xlabel('Time from LEAVE CENTER (s)');
ylabel('LFP (mV)');
title('Trial-averaged LFP (example electrodes)');


%3. now plot the trial-average LFP vs. time split across the reach
%directions. For visualization, just plot the even directions. 
%note: you can (and should) repurpose the code you wrote to make these
%types of figures for spike-rates. 
%FILLIN
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
even_dirs = direction_list(mod(direction_list,2)==0);
plot_electrodes3 = [1 5];
for e=plot_electrodes3
    figure;
    hold on
 for iD=1:length(even_dirs)
     d=even_dirs(iD);
     tr_idx=(trial_reach_direction==d);
     plot(trial_time_lfp, squeeze(mean(trial_lfp(tr_idx,:,e),1)));
 end
 plot([0 0], [min(ylim) max(ylim)], 'k--');
 xlabel('Time from LEAVE CENTER (s)');
 ylabel('LFP (mV)');
 title(['Trial-avg LFP split by reach target (even only), Elec ', num2str(e)]);
end %FILLIN


figs = findall(0,'Type','figure');

for i = 1:length(figs)
    fname = fullfile(dataDir, sprintf('figure_%02d.png', i));
    exportgraphics(figs(i), fname, 'Resolution', 300);
end
