function [cat_hat] = runLeaveOneOutClassification(data, categories)

%[cat_hat] = runLeaveOneOutClassification(data, categories)
%
%function to run leave-one-out cross-validated classification. 
%
%inputs: data - matrix (#trials x #features) of data to be classified
%        categories - vector (#trials x 1) of data category labels
%outputs: cat_hat - vector (#trials x 1) of estimated category for each
%                   trial

num_trials = size(data,1);

cat_hat = nan(num_trials,1); %initialize

%loop through trials
%FILLIN