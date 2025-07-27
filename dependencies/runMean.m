function currAvg = runMean(currAvg, newData, k)
% function for running average. currAvg is the current average and newData
% is the new data that should be added. k is the number of samples in the
% current average + 1. Ignores NaN entries (similar to 'nanmean')

if k <= 0
    error('k has to be at least 1.');
else
    k = round(k); %make sure k is an integer
end

if k == 1
    currAvg = newData; %new average
else
    
    cIdx = ~isnan(currAvg(:)) & ~isnan(newData(:)); %only non-NaN entries
    currAvg(cIdx) = currAvg(cIdx) + ((newData(cIdx) - currAvg(cIdx)) / k); %running average
    
    cIdx = isnan(currAvg(:)) & ~isnan(newData(:)); %new data is non-NaN
    currAvg(cIdx) = newData(cIdx); %add new data to average
    
end
end



