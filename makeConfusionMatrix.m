function confusion_matrix = makeConfusionMatrix(cat_true, cat_hat)

% confusion_matrix = makeConfusionMatrix(class_true, class_hat)
%
%creates a confusion-matrix of the true vs. predicted classified states.
%
%inputs: class_true - vector (#predictions x 1) of true category
%        cat_hat    - vector (# predictions x 1) of estimated category
%
%outputs: confusion_matrix - matrix (#categories x #categories) with probability
%                            of predicting catN when true cat is M for all
%                            N,M pairs
%

category_list = unique(cat_true);      %list of all possible categories
num_categories = length(category_list); %# of categories

confusion_matrix = zeros(num_categories); %initialize confusion matrix

%create confusion matrix
%hint: can loop through each sample and count where each samples goes in the matrix.
%      Or could loop through M and N and count how often that occurred. 
%
%reminder: need to keep track of how often cat_true = M occurs to normalize
%the counts to probabilities. 

%FILLIN


