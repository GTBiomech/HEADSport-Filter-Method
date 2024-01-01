%% HEADSport Filter Method 
% -------------------------------------------------------------------------
% Version: 1.0 (Released: Jan 2024)
% -------------------------------------------------------------------------
% This script calculates the linear (Fh_l) and angular (Fh_a) frequency class 
% using an example case from: 
% Tierney et al. BMJ Open SEM, 2024. http://dx.doi.org/10.1136/bmjsem-2023-001758
% -------------------------------------------------------------------------
% Research contact: Dr Gregory Tierney (g.tierney@ulster.ac.uk)
% -------------------------------------------------------------------------

% Load case
% Note: raw data has already been converted to g and rad/s units 
% GX, GY, GZ represent gyroscope time series data (rad/s); 
% AX, AY, AZ represent accelerometer time series data (g); 
clc; clear; close all;
load('RawData.mat'); 

% iMG Parameters
fs = 3200; % Sampling frequency
range = [200 34.9]; % accelerometer (g) and gyroscope (rad/s) sensor range
T = 1/fs; % Period
L = length(RawSensorData{1,1}.GX); % Data length
Time = (0:L-1)*T'; % Time vector (s)
trigger_time = 0.01; % iMG trigger timepoint (s)                 

% Constraints
Fmax = 320; % Maximum frequency obtained from laboratory head impacts
cfc_min = 60; % Minimum frequency class used in HEADSport_Fh.m
lin_lim = range(1) * 0.95; % 95% of iMG accelerometer range (g)
ang_lim = range(2) * 0.95; % 95% of iMG gyroscope range (rad/s) 

% HEADSport Filter Method
for i = 1:length(RawSensorData)

% Create matrix for linear acceleration (LA) and angular velocity (AV) data
matrix_LA=[RawSensorData{1,i}.AX(1:end),RawSensorData{1,i}.AY(1:end),RawSensorData{1,i}.AZ(1:end)]; 
matrix_AV=[RawSensorData{1,i}.GX,RawSensorData{1,i}.GY,RawSensorData{1,i}.GZ]; 

% Create matrix for artificial (a) raw AV and angular acceleration (AA)
matrix_AVa = CFC(matrix_AV,1/fs,Fmax); 
matrix_AAa = zeros(size(matrix_AVa));
for j = 1:3 
matrix_AAa(:,j) = cent_diff_n(matrix_AVa(:,j),1/fs,5); 
end

% Check iMG range not exceeded based on maximum LA or AV value
lin_max = max(max(abs(matrix_LA))); 
ang_max = max(max(abs(matrix_AV)));
if lin_max > lin_lim || ang_max > ang_lim
    Fh_a = 30; Fh_l =30; 
else
    
% Find column (i.e. X,Y or Z component) associated with max LA or AA value
[~, col] = find(ismember(abs(matrix_LA), max(abs(matrix_LA(:)))));col=col(1); 
[~, col_ang] = find(ismember(abs(matrix_AAa), max(abs(matrix_AAa(:)))));col_ang=col_ang(1); 

% Find pulse duration for LA and AA signal 
[peak_ind, row_pre, row_post] = HEADSport_pulse(matrix_LA, col, trigger_time, fs); 
[peak_ind_ang, row_pre_ang, row_post_ang] = HEADSport_pulse(matrix_AAa, col_ang, trigger_time, fs); 

% Find Power Spectral Density (PSD) characteristics 
[PSD, f, cumsumpsd, LM] = HEADSport_PSD(matrix_LA, row_pre,row_post,col,fs);
[PSD_ang, fa, cumsumpsd_ang, LM_ang] = HEADSport_PSD(matrix_AAa, row_pre_ang,row_post_ang,col_ang,fs);
    
% Find linear (Fh_l) and angular (Fh_a) frequency class
[Fh_l] = HEADSport_Fh(f, cumsumpsd, LM, Fmax, cfc_min);
[Fh_a] = HEADSport_Fh(fa, cumsumpsd_ang, LM_ang, Fmax, cfc_min);
end
end




