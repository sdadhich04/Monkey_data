function [X, Y] = LoadAndPreProcData_filledIn(kinematic_file, spike_file, bin_size)


%function [X, Y] = LoadAndPreProcData(kinematic_file, spike_file, bin_size);
%
%Written by: [YOUR NAME HERE], updated: [DATE HERE]
%
%Loads kinematic and spike data from the specified files and pre-processes
%them for analysis for continuous predictive filters. 
%kinematic data is decimated from the original sampling rate to the
%specified rate (1/bin_size). Spike data is converted from spike-times to
%spike-rates, estimated with bins of size bin_size. 
%
%inputs: kinematic_file - string with full path of kinematic file to load. 
%                         expects this file to contain hand_kinematics and
%                         FS_kinematic variables. 
%         spike_file - string with full path of spike file to load. 
%                      expects this file to contain spike_times cell (#units x 1)
%         bin_size    - double with bin_size (in seconds) to use for firing-rate
%                       estimation and kinematic decimation.
%oututs:  X - kinematic feature matrix [#time bins x #kinematic variables]
%         Y - neural feature matrix [#time bins x #units]
%


%FILLIN
