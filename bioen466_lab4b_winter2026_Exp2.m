%%%%% This script provides an outline for the analyses for lab4b, experiment
%%%%% 2 where we develop and perform initial tests on a Wiener filter.

%%%%% Written by A.L. Orsborn, v200129, updated 211228
%%%%%
%%%%%
%%%%% All lines where you have to fill in information is tagged with a comment including "FILLIN". Use this flag to find everything you need to modify.
%%%%% all figures that need to be included in comprehension questions are
%%%%% flagged with %INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS

%% cell 1: data metadata

% Define some basic things to make it easy to find data files.
% We will want to take advantage of systematic naming structure in our data files.

dataDir = ''; %FILLIN with the path to where your data is stored

codeDir = ''; %FILLIN with the path to where your provided functions for this lab are stored

%we now add all the functions in our code directory to the matlab path
cd(codeDir)
addpath(genpath(codeDir))


%now define which specific file you will analyze
fileBase = 'jeev050712a_bioe466_winter2020';

spikes_file_tag = '_spikes.mat';
kinematics_file_tag = '_kinematics.mat';
lfp_file_tag    = '_lfp.mat';


%% cell 2: Data loading and pre-processing. 
%we will build a winer filter to predict hand velocity (decoded variable) = X
%using unit firing rates (neural feature) = Y. 
%But first, we need to load and prepare the data for this analysis. 

%bin size to use for the filter. 
W_BIN  = 0.2; %bin-size in seconds. 
FS_W   = 1/W_BIN; %sampling rate of wf (HZ = 1/s). 


%we first want to load in our kinematic and neural data and pre-process it for building/testing filters. 
%The pre-processing steps will be to:
%    1) turn our spike counts into rates
%    2) 'bin' our kinematics to have the same time-axis as our binned neural firing rates. 

%1. load the kinematic file and extract the hand velocity
%we will want to match the time-scales of our neural feature and kinematics
%so load the full sampling-rate kinematics. We will down-sample it to match our bin width for firing rates. 

kinematic_file = [dataDir fileBase kinematics_file_tag];
load(kinematic_file, 'hand_kinematics', 'FS_kinematics')  %this notation only loads the variables 'hand_kinematics' and 'FS_kinematics' from the file. 


%to create X, match the sampling rate of X to the bin-size of the data. 
%There are many ways to down-sample data. We will use the matlab function
%'decimate', which approximates the data via a low-pass filter (rather than
%just taking very Nth data point, which can be noisy. 
%Look at the help for that function and understand how to use it to create
%an X vector sampled at 1/W_BIN rate. Hint: think about the units of
%FS_kinematics and W_BIN

num_kin = size(hand_kinematics,2); %there are multiple kinematic variables. 

downsample_ratio = floor(); %FILLIN. We use 'floor' to make the ratio an integer number. Fill in the argument to pass into floor. 

X = nan(ceil(size(hand_kinematics,1)/downsample_ratio), num_kin); %initialize X [#time-points of down-sampled data by num_kin]
%note: based on the help of 'decimate', we know exactly what size to expect for X

%loop through each kinematic variable (columns of hand_kinematics) to decimate each one (decimate does not work on matrices)
for i=1:num_kin
    X(:,i) = decimate(); %FILLIN
end

%create a time-axis for the full and down-sampled data so we can visualize them.  
time_kin = (0:size(,1)-1)./FS_kinematics; %FILLIN
time_x = (0:size(,1)-1)./FS_W; %FILLIN


%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure
plot(time_kin, hand_kinematics(:,1)) %plot just one of the variables to see what's happening
hold on
plot(time_x, X(:,1), '--') %make the filtered versions a dashed line 
xlabel('')
ylabel('')
set(gca, 'xlim', [10 20]) %zoom in some so you can see

%Look at the plot and make sure you understand what this operation has done. 
%Do you notice anything about the temporal relationship between the
%decimated and non-decimated versions of your data? Where does that come
%from?


%2. load the spike file and estimate your firing rates using the bin-size
%of your winer filter. 

spike_file = [dataDir fileBase spikes_file_tag];
load(spike_file, 'spike_times') %this notation only loads the variable 'spike_times' from the file. 

num_units = length(spike_times);

%now, we want to bin our spikes. We will be running continuous predictions
%(i.e. not aligning our data in time based on task events). 
%our past binning code, however, operated on trial-aligned data. 
%But! We can repurpose this code to still achieve our goals. 
%We can consider our stream of spike-times as coming from a single trial
%that is the length of the file. Our 'binTrialAlignedSpikes.m' (lab 2) will 
%still work. 

%create a 'spike_labels' variable (to use with our bin function)
%for each unit to specify that all
%spikes belong to the same trial (trial #1)
spike_labels = cell(size(spike_times));
for iU=1:num_units
    spike_labels{iU} = ones(size()); %FILLIN
end

%now bin our spikes. 
%Since we're treating the file as one continuous trial:
%    time_before = 0
%    time_after = total time of the file
%We can use the hand_kinematics and FS_kinematics to find the total file duration.
%
time_before = 0;
time_after  = size(,1)/(); %FILLIN

[spike_rate, time_bins] = binTrialAlignedSpikes_filledIn(...
    spike_times, spike_labels, time_before, time_after, W_BIN);

%we can get rid of the 1st dimension of spike_rate because it's 1 (for the single trial)
spike_rate = squeeze(spike_rate); %spike_rate now [time x #units] vector


%Now set Y for our wiener filter
Y = ; %FILLIN. hint: it's a variable you already have. 

%double-check that X and Y have the same time-axis. 
%the histogramming and decimation procedures we do to bin data often leaves us with 1
%less time-point for one of the matrices (depending on rounding done with decimation).
%We will remove any extra data point(s) and press on. 
time_length = min( size(Y,1), size(X,1) );
diff_y = size(Y,1)-time_length;
diff_x = size(X,1) - time_length; 
if diff_y > 0
    disp(['Trimming ' num2str(diff_y) ' time-points from Y'])
    %note: if diff_y > 1, something is probably off with your code. 
end
if diff_x > 0
    disp(['Trimming ' num2str(diff_x) ' time-points from X'])
    %note: if diff_x > 1, something is probably off with your code.
end
Y = Y(1:time_length,:);
X = X(1:time_length,:);



%These operations (loading the kinematics and spikes; pre-processing them with binning)
%will be used repeatedly in our analysis for this lab. 
%To clean up our code, let's make a function that does these operations.  
%
%The inputs to the function will be the kinematic and spike files to load, and the bin-size. 
%The outputs will be the processed kinematic (X) and spike data (Y) 
%
%we've provided a shell function LoadAndPreProcData.m. Open it and fill in. 
%Important: You can COMPLETELY repurpose the code you wrote above. Just turn it into a
%function. (You won't need to plot the data in the function--that's just
%for our understanding.). 


%% cell 3: Now we want to actually train our filter.
%we will begin with a 1-lag filter (i.e. using a single time-point) to understand the steps, and then
%use a more general function for training a filter with multiple lags. 

%clean up some variables
clear spike_* diff_* time_* X Y

%specify parameters of the filter (bin size and # lags.) 
W_BIN  = 0.2; %bin-size in seconds. 
W_LAGS = 1;

%data
spike_file = [dataDir fileBase spikes_file_tag];
kinematic_file = [dataDir fileBase kinematics_file_tag];


%use our preprocessing function to get the neural and kinematic data. 
[X, Y] = LoadAndPreProcData(); %FILLIN

%We will only predict the hand velocity (x,y), which is a subset of the
%kinematic feature matrix X. Sub-select the columns of data that correspond
%to the hand velocity:
X = X(:,[]); %FILLIN


%We will also add a constant state to Y. This allows the states to have a
%non-zero mean. This corresponds to the b term in X(t) = b + sum( a_i*Y_i(t)) . 
Y(:,end+1) = ones(size(Y,1),1); %this 'end+1' notation is a matlab short-hand to add new data. This adds a column of ones to the end of X. 


%Now, let's train a Wiener filter. 

%The wiener model is that Y = A*X
%The optimal (in a least-squares sense) A can be learned from training data
%X and Y:
%A = inv(Y'*Y)*Y'*X;
%note: in matlab a more efficient way to invert matrices is with the '\' operator. 
%x = inv(a)*b is roughly equivalent to x = a\b. However, they're computed differently. 
%Matlab recommends the \ method. Use that here. 
A = ; %FILLIN

%Now use your learned A matrix to generate your predictions, X_hat:
X_hat = ; %FILLIN


%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
%plot X vs. X-hat for all your variables (X and Y)
figure
y_axis_labels = {'X velocity (cm/s)', 'Y velocity (cm/s)'};
for i=1:2
    subplot(2,1,i)
    plot(X(:,i))
    hold on
    plot(X_hat(:,i), 'r')
    xlabel()
    ylabel(y_axis_labels{i})
    set(gca, 'xlim', [1 60/W_BIN]) %zoom-in to plot 1 minute (60sec) of data
    title('Predictions on training data')
end

%finally, let's compute the coefficient of determination for our data. 
%we've provided a function coeffDetermination to perform these
%calculations. Look at the help and use it to calculate Rsquared. 
Rsquared = coeffDetermination(); %FILLIN

disp(['R^2 for X velocity = ', num2str(Rsquared(1))])
disp(['R^2 for Y velocity = ', num2str(Rsquared(2))])



