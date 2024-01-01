function [Fh] = HEADSport_Fh(f, cumsumpsd, LM, Fmax, cfc_min)
% -------------------------------------------------------------------------
% Version: 1.0 (Released: Jan 2024)
% -------------------------------------------------------------------------
% This script calculates the frequency class (Fh) based on HEADSport filter method, see Figure 3 from: 
% Tierney et al. BMJ Open SEM, 2024. http://dx.doi.org/10.1136/bmjsem-2023-001758
% -------------------------------------------------------------------------
% Research contact: Dr Gregory Tierney (g.tierney@ulster.ac.uk)
% -------------------------------------------------------------------------
% Input variable:
% f: Frequency lines
% cumsumpsd: Normalised (%) cumulative sum of PSD for frequency lines
% LM: Index of local minima frequencies from PSD
% Fmax: Max frequency from laboratory impacts
% cfc_min: Min frequency class used in HEADSport_Fh.m
% -------------------------------------------------------------------------
% Output variable: 
% Fh = Frequency Class
% -------------------------------------------------------------------------

    % Step 1: Get percentile PSD
    ind1=find(cumsumpsd >=0.95,1,'first');% Index of frequecny at which 95% PSD (cumulative sum) is achieved
    Fh = round(f(ind1)); % Round frequency to nearest interger
    if Fh == 0; Fh = 1000;end % Make Fh exceed Fmax (see Step 3)
    if Fh<cfc_min; Fh=cfc_min; end % Fh is cfc_min
   
    % Step 2: Check first local minima frequency if Fh<Fmax. Fh is lowest of the
    % two values
    if Fh <= Fmax 
        try
        indmi = find(LM == 1); % Find index of first local minima
        fh_mi = round(f(indmi(1)));  % Frequency associated with minima (rounded to nearest integer)
        if fh_mi < Fh; Fh = fh_mi; end % Fh is lowest of the two values
        if Fh < cfc_min; Fh = cfc_min; end % Fh is cfc_min
        end
    end
    
    % Step 3: Check first local minima frequency if Fh>Fmax. Use it if < Fmax 
    if Fh > Fmax 
        try
        ind1 = find(LM == 1); % Find index of first local minima
        Fh = round(f(ind1(1))); % Round frequency to nearest interger
        if Fh < Fmax; end % Fh is first local minima
        if Fh < cfc_min; Fh = cfc_min; end % Fh is cfc_min
        catch % No minima
        Fh = cfc_min; % Fh is cfc_min
        end
    end
