%%%%% This script provides an outline for the analyses for lab4, experiment
%%%%% 2 where we develop and test classifiers to predict movement direction
%%%%% from neural activity.

%%%%% Written by A.L. Orsborn, v191227
%%%%%
%%%%%
%%%%% All lines where you have to fill in information is tagged with a comment including "FILLIN". Use this flag to find everything you need to modify.
%%%%% all figures that need to be included in comprehension questions are
%%%%% flagged with %INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS

%% data metadata

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


%define some constants that help with data interpretation
GOCUE_CODE      = 5; %event-code flagging when the go-cue comes on
LEAVE_CENTER_CODE= 6; %event-code flagging when hand leaves center
REWARD_CODE     = 9; %event-code flagging when reward turns on
AT_TARGET_CODE    = 7; %event-code flagging when the animal is at the target
TIMEOUT_CODE      = 8; %event-code flagging when a reach time-out (didn't get to target fast enough) error occurred

TARGET_OFFSET   = 63; %offset to subtract from event-codes for target direction


%%

%%%%%%%% A: Building a classifier and visualizing/quantifying it's
%%%%%%%% performance.

% We will first pick a single data set (neural feature, amount of training
% data, etc.) and train a linear-discriminant classifier.
%Later portions of this lab will explore how our decoding performance
%varies as we change these properties.


% We will start with spike firing rates as our 'neural feature'.
% 1. Get trial-aligned neuron firing rates
% based on your exploration of the data in experment 1, pick the alignment-event
% that you think would give the best decoding performance (hint: you want
% to minimize trial-to-trial variability)
% load enough data before/after to explore how well decoding varies over
% time. (2 seconds before/after should be sufficient)

spike_file = ; %FILLIN
time_before = ; %FILLIN. in seconds
time_after = ;  %FILLIN. in seconds
bin_width  = ;   %FILLIN. in seconds
align_code = ;   %FILLIN



%2. get the "labels" for each trial that we will decode.
%   We will decode the reach target identity (reach_direction)
%   reminder: event-information is saved in spike_file.
%FILLIN


%3. Pick a time-window to use for decoding.
%   neural 'feature' is the mean firing in this time-window.
predict_window = [-.4 .4]; 
time_inds = ; %FILLIN. A logical vector that is 1 for all times within your prediction window and 0 otherwise. 

neural_feature = ; %FILLIN. get average over selected time-window. Should be trials x units


%4. Divide the data into a training set and a test-set
num_trials = size(neural_feature,1);
fraction_data_train = 0.5;

%create a logical that is true when data should be used for training, false
%otherwise. Randomly select the trials to use/leave-out
use_for_training = false(num_trials,1);
rand_trial_order = randperm(num_trials);
use_for_training(rand_trial_order(1:round(num_trials*fraction_data_train))) = true;

%5. use the 'classify' function in matlab to predict reach_direction_hat
%for the left-out data (i.e. where use_for_training = 0)
%We will use a linear classifier.
%FILLIN


%6. Compute, overall, the probability of correct prediction for your classifier.
p_correct = ; %FILLIN

fprintf('Overall prediction accuracy = %4f\n', p_correct)


%7. create a confusion matrix and plot it.
%   we have provided a shell function makeConfusionMatrix.m for you to fill in
confusion_matrix = makeConfusionMatrix(reach_direction_true, reach_direction_hat);

%visualize the confusion matrix
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure
imagesc(confusion_matrix, [0 1])
xlabel('Predicted Target')
ylabel('True target')
colormap hot
colorbar
title('Confusion matrix of reach target decoding')


% Try re-running this analysis again to see how it varies with the sub-set of data used for training/testing
% What technique we discussed allows you to estimate this variance?
% Also try re-running this varying the fraction of data used for
% training to understand how the amount of training data influences
% results.


%% Now let's do a slight modification to our decoding to do leave-one-out cross-validation.
%  This technique uses all trials except one to train, and then tests on
%  the left-out trial. To get an overall accuracy across the entire
%  dataset, we do this procedure for every trial.

%We've provided a shell function for you to fill-in to do leave-one-out
%cross-validation
reach_direction_hat = runLeaveOneOutClassification(neural_feature, reach_direction);

%calculate your overall accuracy and plot a confusion matrix.

p_correct = ; %FILLIN

fprintf('Overall prediction accuracy = %4f\n', p_correct)

confusion_matrix = makeConfusionMatrix(reach_direction, reach_direction_hat);

%visualize the confusion matrix
figure
imagesc(confusion_matrix, [0 1])
xlabel('Predicted Target')
ylabel('True target')
colormap hot
colorbar
title('Confusion matrix of reach target decoding')


%%
%%%%%%%% B: test how classifier performance depends on the amount of neural
%%%%%%%% data used.

%now that we have the basic operations down, let's explore how decoding
%varies depending on the data we put in.


%1. Generate a plot of how prediction accuracy varies as a function of the
%number of units used. Use leave-one-out cross-validation.
%When you use less than the total number of units, remove units randomly.
%on your plot, draw a horizontal line to indicate the maximum decoding
%performance and a vertical line at the window-length it corresponds to
%FILLIN
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS

%%
%2. Generate a plot of how prediction accuracy varies as a function of the
%predict_window used. Use leave-one-out cross-validation.
%Keep the window START fixed at -0.5 and vary the total length of the
%window from [0.1 to 1.5] seconds.
%on your plot, draw a horizontal line to indicate the maximum decoding
%performance and a vertical line at the window-length it corresponds to

%FILLIN
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
%%

%Now let's explore how the neural features themselves influence decoding.

%3. Generate a plot of how prediction accuracy varies across each individual unit.
% pick a predict_window that seems to work well from your previous results.

%FILLIN



%%

%4. Generate a plot (heatmap) of how prediction accuracy varies across each individual unit
%   AND the time-window used
% specifically: keep the prediction window LENGTH fixed at 0.5. For each
% unit, sweep window start-times ranging from -1.5 to 1

%FILLIN



%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure
imagesc() %FILLIN
colormap hot
colorbar
ylabel('Unit ID')
xlabel('Prediction Window Start (s)')
title('Prediction accuracy vs. time and unit, window length = 0.5')

%%

%5. Finally, try using a different neural feature other than spike firing rates.
% Let's use the LFP voltages as a neural feature.

% 1. Get trial-aligned LFP voltages
% based on your exploration of the data in experment 1, pick the alignment-event
% that you think would give the best decoding performance (hint: you want
% to minimize trial-to-trial variability)
% load enough data before/after to explore how well decoding varies over
% time. (2 seconds before/after should be sufficient)



time_before = ; %FILLIN. in seconds
segment_length = ;  %FILLIN. in seconds
align_code = ; %FILLIN

align_times = ; %FILLIN
[trial_lfp] = trialAlignData(); %FILLIN
trial_time_lfp = ; %FILLIN


%2. "labels" for each trial are already loaded.


%3. Pick a time-window to use for decoding.
%   neural 'feature' is the mean in this time-window.
predict_window = [-.4 .4]; 
time_inds = ; %FILLIN. Logical vector that is 1 when time is within predict_window and 0 otherwise. 

neural_feature = ; %FILLIN. get average over selected time-window. Should be trials x units


%4. Divide the data into a training set and a test-set
num_trials = size(neural_feature,1);
fraction_data_train = 0.5;

%create a logical that is true when data should be used for training, false
%otherwise. Randomly select the trials to use/leave-out
use_for_training = false(num_trials,1);
rand_trial_order = randperm(num_trials);
use_for_training(rand_trial_order(1:round(num_trials*fraction_data_train))) = true;

%5. use the 'classify' function in matlab to predict reach_direction_hat
%for the left-out data (i.e. where use_for_training = 0)
%We will use a linear classifier.
%FILLIN

%6. Compute, overall, the probability of correct prediction for your classifier.
p_correct = ; %FILLIN

fprintf('Overall prediction accuracy = %4f\n', p_correct)

%7. create a confusion matrix and plot it.
confusion_matrix = makeConfusionMatrix(reach_direction_true, reach_direction_hat);

%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure
imagesc(confusion_matrix, [0 1])
xlabel('Predicted Target')
ylabel('True target')
colormap hot
colorbar
title('Confusion matrix of reach target decoding')


% As time allows, try varying the amount of data (prediction_window, number of channels) for the LFP signals.
%Do you expect the prediction accuracy vs. prediction_window curve to be the same or different for this neural feature? Why?
%Do you expect the prediction accuracy vs. number of channels curve to be the same or different for this neural feature? Why?
%What do differences in these curves tell you about how reach target information is "encoded" in each neural feature?
