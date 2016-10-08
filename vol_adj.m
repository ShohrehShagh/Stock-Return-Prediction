function [X_vol_adj, vol] = vol_adj(X, com1, com2, min_periods, vol_delay)

% Volatility estimation and adjustment
% This function takes in a matrix X along with EWMA parameters and returns
% the estimated volatilities for columns of X as well as the
% volatility-adjusted values of the columns of X

% INPUTS:
%   X          : The input matrix with columns corresponding to time series
%   com1        : center of mass parameter for variance EWMA
%   com2        : centre of mass parameter for mean EWMA
%   min_periods : minimum number of non-nan values before returning a 
%                non-nan EWMA value, and volatility
%   vol_delay  : the delay in using the estimated volatilities for
%                adjustment (should be at least 1)
% Either of these parameters, except X, can be set to [], in which case, 
% they will be set to their default values.

% OUTPUTS:
%   vol       : Estimated volatilities
%   X_vol_adj : volatility-adjusted values

% Setting the parameters not provided to their default values
if nargin<2
    com1 = []; com2 = []; min_periods = []; vol_delay = [];
end;
if nargin==2
    com2 = []; min_periods = []; vol_delay = [];
end;
    
    
if isempty(com1)
    com1 = 20;
end
if isempty(com2)
    com2 = 120;
end

if isempty(min_periods)
        min_periods = 3*com1;
end
if isempty(vol_delay)
    vol_delay = 1;
end

vol = NaN(size(X));
X_vol_adj = NaN(size(X));

for ii=1:size(X,2)
    cur = X(:,ii);

    if ~isnan(com2)
        mncur = EWMA(cur,'com',com2,'min_periods',min_periods);
        cur = cur - mncur;
    end
    vol(:,ii) = sqrt(EWMA(cur.^2,'com',com1,'min_periods',min_periods));
    X_vol_adj(:,ii) = [NaN(vol_delay,1); X((vol_delay+1):end,ii)./...
        vol(1:(end-vol_delay),ii)];
end
