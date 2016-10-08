function y = EWMA(x,varargin)

% Exponentially Weighted Moving Average calculation
% Based on Python's 'ewma' function

% Decaying weights can be defined in a few ways:
%   centre of mass (com): \alpha = 1/(1+com)
%   halflife: \alpha = 1-exp(log(0.5)/halflife) 
%   span: \alpha = 2/(span+1)

% This implementation does not support resampling of the data

% Inputs
% Required: 
% x        : vector to be averaged
% Optional:  
% com         : float
%               Specify decay in terms of center of mass
% span        : float
%               Specify decay in terms of span
% halflife    : float, optional
%               Specify decay in terms of halflife
% min_periods : int, default 0
%               Minimum number of observations in window required to have
%               a value (otherwise result is NA).
% adjust      : boolean, default True
%               Divide by decaying adjustment factor in beginning periods
%               to account for imbalance in relative weightings
%               (viewing EWMA as a moving average)
% ignore_na   : boolean, default False
%               Ignore missing values when calculating weights
%

p = inputParser;

% Setting default values of parameters
defaultCom = NaN;
defaultSpan = NaN;
defaultHalflife = NaN;
defaultAdjust = true;
defaultMin_periods = 0;
defaultIgnore_na = true;

% Adding parameter name-value pairs argument to input parser scheme
addRequired(p,'x',@isnumeric);
addParameter(p,'com',defaultCom,@isnumeric);
addParameter(p,'span',defaultSpan,@isnumeric);
addParameter(p,'halflife',defaultHalflife,@isnumeric);
addParameter(p,'adjust',defaultAdjust,@islogical);
addParameter(p,'min_periods',defaultMin_periods,@isnumeric);
addParameter(p,'ignore_na',defaultIgnore_na,@islogical);

% Parsing function inputs
parse(p,x,varargin{:});
min_periods = p.Results.min_periods;
adjust = p.Results.adjust;
ignore_na = p.Results.ignore_na;

x = x(:);
N = length(x);

% Calculate the center of mass (transforming span and halflife to com)
com = get_centre_of_mass(p.Results.com,p.Results.span,p.Results.halflife);

% Calculating the decay parameter from center of mass
alpha = 1/(1 + com);

minp = max(min_periods,1);
old_wt_factor = 1 - alpha;

if (adjust) 
    new_wt = 1; 
else
    new_wt = alpha; 
end

weighted_avg = x(1);
y = NaN(size(x));

is_observation = ~isnan(x);
nobs = cumsum(is_observation);

if nobs(1) >= minp
    y(1) = weighted_avg; 
end
old_wt = 1;    

for ii = 2:N
    cur = x(ii);           
    if ~isnan(weighted_avg)        
        if (is_observation(ii)) || (~ignore_na) 
            old_wt = old_wt * old_wt_factor; % updating the weight
            if (is_observation(ii))
                if weighted_avg ~= cur
                    % updating the weighted average:
                    weighted_avg = ((old_wt * weighted_avg) + ...
                        (new_wt * cur)) / (old_wt + new_wt); 
                end
                if adjust
                    old_wt = old_wt + new_wt;
                else 
                    old_wt = 1;
                end
            end
        end    
    elseif is_observation(ii)
        % If no value is assigned to 'weighted_avg', set it to 'cur'
        weighted_avg = cur;
    end
    if nobs(ii) > minp
        % output 'weighted_avg' if we have seen enough non-nan values in x
        y(ii) = weighted_avg;
    end
end    

function com_out = get_centre_of_mass(com,span,halflife)
% This function transforms the given decay parameter to center of mass

numvalid = ~isnan(com)+ ~isnan(span) + ~isnan(halflife);
if numvalid < 1
    % No decay parameter specified
    msgID = 'ewma: no weight';
    msg = 'No weight specified for ewma.';
    baseException = MException(msgID,msg);
    throw(baseException);
end
if numvalid > 1
    % More than one weight specified
    msgID = 'ewma: multiple weights';
    msg = 'Multiple weights specified for ewma.';
    baseException = MException(msgID,msg);
    throw(baseException);
end

% Calculating center of mass
if ~isnan(span)
    com_out = (span-1)/2;
elseif ~isnan(halflife)
    decay = 1- exp(log(0.5)/halflife);
    com_out = 1/decay - 1;
else
    com_out = com;
end    


