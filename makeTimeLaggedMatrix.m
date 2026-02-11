function Y2 = makeTimeLaggedMatrix(Y, num_lags)

%Y2 = makeTimeLaggedMatrix(Y, num_lags)
%
%creates a time-lagged matrix of data for training/prediction with a multi-lag Wiener
%filter. Also appends a column of 1's, which allows modeling of a DC offset in predictions. 
%
%inputs: Y - [#time-points x #measurements] matrix 
%        num_lags - #lags to create
%
%outputs: Y2 - [#time-points - num_lags+1 x #measurements*num_lags+1]
%matrix where time-shifted versions of the measurement data have been
%appended to columns.

[Nt, Nspk] = size(Y);

%Initialize (note that this appends a column of ones to the end)
Y2 = ones(Nt-num_lags+1, Nspk*num_lags + 1);
 
%loop through lags
for iL=1:num_lags
    Y2(:, (iL-1)*Nspk+1:iL*Nspk) = Y( (1:Nt-num_lags+1)+iL-1, :);
end