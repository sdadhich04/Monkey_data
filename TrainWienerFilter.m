function A = TrainWienerFilter_filledIn(X, Y, num_lags)


%A = trainWiener(X, Y, num_lags)    (A. Orsborn, created 6-16-10, updated 1/20/20)
%
%Trains wiener filter to predict states (X) based on observations (Y)
%Uses ALL input data in matrices X and Y for training.
%
%inputs:  Y     - observations i.e. binned neural activity (time x #observations)
%         X     - state matrix i.e. binned kinematics      (time x #states)
%                 will fit matrices for all input states
%         num_lags- # time lags to use
%outputs: A     - trained matrix (#observations*num_lags x #states)
%
 
 
%get time length & # observations (i.e. # neurons)
[N_time, N_obs] = size(Y);
 
%double-check for consistency between X and Y, # lags
if size(X,1) ~= N_time
    error('X and Y matrices must have same # of rows (i.e. time).\nsize(X) = %g  %g\nsize(Y) = %g %g', ...
        N_time, N_obs, size(X,1), size(X,2))
end
if num_lags>N_time
    error('num_lags > time length!')
end
 
 
 
%%%%%%%%training %%%%%%%
 
%tile Y to create time x (neurons*lags) training matrix
%hint: use the provided function makeTimeLaggedMatrix
Y2 = makeTimeLaggedMatrix(Y, num_lags);

 
 
%calculate A matrix ( A = inv(Y'*Y)*Y'*X )
A = ; %fill-in
  
end