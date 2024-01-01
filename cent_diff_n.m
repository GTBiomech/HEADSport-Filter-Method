function df = cent_diff_n(f,h,n)

% Cite as: James (2023). N-Point Central Differencing (https://www.mathworks.com/matlabcentral/fileexchange/36123-n-point-central-differencing), MATLAB Central File Exchange.
% df = cent_diff_n(f,h,n)
% Computes an n-point central difference of function f with spacing h.
% Returns a vector df of same size as f.
% Input f must be a vector with evenly spaced points.
% Input n must be 3,5,7, or 9.
% All three inputs are required.
%
% Differences for points near the edges are calculated with lower order.
% For example, if n=5 and length(f)=10, then 3-point central differencing is used
% to calculate values at points 2 and 9, 2-point forward differencing is used for
% point 1, 2-point backward differencing is used for point 10, and 5-point central
% differencing is used for points 3-7.
%
% If f contains less than n points, the order will be downgraded to the
% maximum possible.  Ex: if length(f) = 6, n will be downgraded to 5.
%
% Differencing formulae from: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/


if nargin < 3
    error('Not enough inputs.  See help documentation.')
end

if ~isscalar(h)
    error('Input h must be a scalar value.')
end

    
possible_ns = [3,5,7,9];
if ~ismember(n,possible_ns)
    error('Input n must be 3,5,7, or 9.')
end

numPts = length(f);
if numPts < n
    newN = max(possible_ns(possible_ns<=numPts));
    warnstr = [num2str(n) '-point differencing was requested,\n'...
       'but input function only has ' num2str(numPts) ' points.\n'...
       'Switching to ' num2str(newN) '-point differencing.'];
    warning(warnstr,'%s')
    n = newN;
end

df_1 = b_diff(f,h);
df_End = f_diff(f,h);

% Calculate 3-point for all
df_3pt = c_diff(f,h,3);

if n >=5
    df_5pt = c_diff(f,h,n);
    % For the 2nd and next-to-last grid point, use 3-point differencing.
    df_2 = df_3pt(1); 
    df_Endm1 = df_3pt(end);
end
if n >=7
    df_7pt = c_diff(f,h,7);
    % For the 3nd and 2nd from last grid point, use 5-point differencing.
    df_3 = df_5pt(1);
    df_Endm2 = df_5pt(end);
end
if n>= 9
    df_9pt = c_diff(f,h,9);
    % For the 4nd and 3rd from last grid point, use 7-point differencing.
    df_4 = df_7pt(1);
    df_Endm3 = df_7pt(end);
end

switch n
    case 3 
        df = [df_1 df_3pt df_End];
    case 5
        df = [df_1 df_2 df_5pt df_Endm1 df_End];
    case 7
        df = [df_1 df_2 df_3 df_7pt df_Endm2 df_Endm1 df_End];
    case 9
        df = [df_1 df_2 df_3 df_4 df_9pt df_Endm3 df_Endm2 df_Endm1 df_End];
end

        


end

function df = c_diff(f,h,n)
midStartPoint = ceil(n/2); % First point at which full n points can be used
midEndPoint = length(f)-midStartPoint+1; % Last point at which full n points can be used

df = [];
for k = midStartPoint:midEndPoint
    switch n
        case 3
            df_k = (f(k+1) - f(k-1))/(2*h);
        case 5
            df_k = (f(k-2) - 8*f(k-1) + 8*f(k+1) - f(k+2))/(12*h);
        case 7
            df_k = (-f(k-3) + 9*f(k-2) - 45*f(k-1) + 45*f(k+1) - 9*f(k+2) + f(k+3))/(60*h);
        case 9
            df_k = (3*f(k-4) - 32*f(k-3) + 168*f(k-2) - 672*f(k-1) + 672*f(k+1) - 168*f(k+2) + 32*f(k+3) - 3*f(k+4))/(840*h);
    end
    df = [df df_k];
end
end

function df1 = b_diff(f,h)
df1 = (f(2)-f(1))/h;
end

function dfEnd = f_diff(f,h)
dfEnd = (f(end)-f(end-1))/h;
end