function [peak_ind, row_pre, row_post] = HEADSport_pulse(matrix, col, trigger_time, fs)
% -------------------------------------------------------------------------
% Version: 1.0 (Released: Jan 2024)
% -------------------------------------------------------------------------
% This script calculates kinematic pulse characteristics from: 
% Tierney et al. BMJ Open SEM, 2024. http://dx.doi.org/10.1136/bmjsem-2023-001758
% -------------------------------------------------------------------------
% Research contact: Dr Gregory Tierney (g.tierney@ulster.ac.uk)
% -------------------------------------------------------------------------
% Input variable:
% matrix: Data structure containing kinematics
% col: Column (X,Y or Z) associated with max LA or AA value
% trigger_time: iMG trigger timepoint (s)
% fs: Sampling frequency
% -------------------------------------------------------------------------
% Output variable: 
% peak_ind: Index of peak value
% row_pre: First zero crossing index (zci) before peak
% row_post: First zero crossing index (zci) after peak
% -------------------------------------------------------------------------

[~,peak_ind] = max(abs(matrix(:,col))); % Get index of peak value

% Use zero crossing (zc) to define pulse i.e., where signal passes 0 before and after peak 
v=matrix(:,col); % Column with highest magnitude
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);  % Find zero crossings
zx = zci(v); % Get index

% Find first zero crossing index (zci) before peak
try   
zci_ind_bef=find(zx<peak_ind);row_pre=zx(zci_ind_bef(end));
catch
row_pre=trigger_time*fs - 0.01*fs + 1; % Index associated with trigger time minus 10 ms (requires at least 10 ms of data before trigger)
end

%find first zero crossing index (zci) after peak
try
   zci_ind_aft=find(zx>peak_ind); row_post=zx(zci_ind_aft(1)); 
catch
   row_post=trigger_time*fs + 0.04*fs; % Index associated with trigger time plus 40 ms (requires at least 40 ms of data after trigger)
end