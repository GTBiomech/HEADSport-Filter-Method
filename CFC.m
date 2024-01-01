function CFC_out=CFC(Data,T,cfc)
% -------------------------------------------------------------------------
% Channel Frequency Class Filter
% -------------------------------------------------------------------------
% Description given in:
% Tierney et al. BMJ Open SEM, 2024. http://dx.doi.org/10.1136/bmjsem-2023-001758
% -------------------------------------------------------------------------
% Input variable:
% T: Time Step (s)
% Data: Data structure containing kinematics to be filtered
% cfc: Frequency class
% -------------------------------------------------------------------------
% Output variable: 
% CFC_out: Data structure containing filtered kinematics
% -------------------------------------------------------------------------

wd=2*pi*cfc*2.0775;
wa=(sin(wd*T/2))/(cos(wd*T/2));
a0=(wa^2)/(1+wa*(2^.5)+wa^2);
a1=2*a0;
a2=a0;
b1=-2*((wa^2)-1)/(1+wa*(2^.5)+wa^2);
b2=(-1+wa*(2^.5)-wa^2)/(1+wa*(2^.5)+wa^2);

A=[1;-b1;-b2];
B=[a0;a1;a2];

CFC_out=filtfilt(B,A,Data);
