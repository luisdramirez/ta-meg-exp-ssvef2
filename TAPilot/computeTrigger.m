function trigger = computeTrigger(varargin)

a = cell2mat(varargin);

trigger = nansum(2.^(a-1));

if isnan(trigger)
    trigger = 0;
end