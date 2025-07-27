function out = insertColumns(X, cIdx, cVal, cWidth)
%function to insert columns filled with 'cVal' at the locations defined by
% 'cIdx' into the matrix X. cWidth defines the width of each inserted
% column.
% smusall, 02.2022

if ~exist('cWidth', 'var') || isempty(cWidth)
    cWidth = 1;
end

if islogical(cIdx)
    cIdx = find(cIdx);
else
    cIdx = sort(cIdx, 'ascend');
end
cIdx(cIdx < 1) = []; %ignore non-positive indices
cIdx(cIdx > size(X,2)) = []; %ignore indices that are too large

% make index for where to put data into a larger array that has columns
dataIdx = false(1, size(X,2) + (length(cIdx)*cWidth));
cIdx = [1, cIdx];

for x = 1 : length(cIdx)
    if x == length(cIdx)
        dataIdx(cIdx(x) : end) = true;
    else
        dataIdx(cIdx(x) : cIdx(x+1)-1) = true;
        cIdx = cIdx + cWidth;
    end
end
    
dSize = size(X);
dSize(2) = length(dataIdx); %new size with added columns
out = repmat(cVal, dSize); %initialize output array
out(:, dataIdx) = X;
