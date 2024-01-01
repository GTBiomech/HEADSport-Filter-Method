function [PSD, f, cumsumpsd, LM] = HEADSport_PSD(matrix, row_pre, row_post, col, fs)
% -------------------------------------------------------------------------
% Version: 1.0 (Released: Jan 2024)
% -------------------------------------------------------------------------
% This script calculates the Power Spectral Density (PSD) Characteristics from: 
% Tierney et al. BMJ Open SEM, 2024. http://dx.doi.org/10.1136/bmjsem-2023-001758
% -------------------------------------------------------------------------
% Research contact: Dr Gregory Tierney (g.tierney@ulster.ac.uk)
% -------------------------------------------------------------------------
% Input variable:
% matrix: Data structure containing kinematics
% row_pre: First zero crossing index (zci) before peak
% row_post: First zero crossing index (zci) after peak
% col: Column (X,Y or Z) associated with max LA or AA value
% fs: Sampling frequency
% -------------------------------------------------------------------------
% Output variable: 
% PSD: Power Spectral Density (PSD) associated with frequency lines
% f: Frequency lines
% cumsumpsd: Normalised (%) cumulative sum of PSD for frequency lines
% LM: Index of local minima frequencies from PSD
% -------------------------------------------------------------------------
    
L=length(matrix(row_pre:row_post,col)); % Data length 
NFFT = 2^nextpow2(L); % Finds length of padded signal that is input into fft, see: https://uk.mathworks.com/matlabcentral/answers/440116-how-do-i-use-fft-to-calculate-the-psd-of-a-signal-that-is-padded-with-zeros
xdft  = fft(matrix(row_pre:row_post,col),NFFT); % Returns FFT of signal for the padded input length
xdft = xdft(1:NFFT/2+1); % Single sided
PSD = (1/(fs*L)) * abs(xdft).^2; % Noting that the zeros used to pad the input do not contribute to the total time of the signal and so 'N', the number of points in the original signal, and not 'NFFT' is used 
PSD(2:end-1) = 2*PSD(2:end-1); % Contribution of the values that were truncated in the previous step can be included
f = 0:fs/NFFT:fs/2; % Frequency lines
cumsumpsd = cumsum(PSD)./sum(PSD); % Calculate normalised (%) cumulative sum of PSD for frequency lines
LM = islocalmin(PSD); % Index of local minima frequencies from PSD
