function DataOut = arrayCrop(DataIn,mask)
%code to apply a mask to an image or stack of images. This needs a single
%or double as the datatype for 'DataIn'. The output will converted to
%double otherwise.

dSize = size(DataIn); %size of input matrix
if dSize(1) == 1
    DataIn = squeeze(DataIn); %remove singleton dimensions
    dSize = size(DataIn);
end

if length(dSize) == 2
    if dSize(1) == 1
        DataIn = DataIn';
        dSize = size(DataIn); %size of input matrix
    end
    dSize = [dSize 1];
end
mSize = size(mask);


%check if datatype is single. If not will use double as a default.
if isa(DataIn,'single')
    dType = 'single';
else
    dType = 'double';
end


if dSize(1) > size(mask,1) || dSize(2) > size(mask,2)
    DataIn = DataIn(1:size(mask,1), 1:size(mask,2), :); %cut to size
end
if dSize(1) < size(mask,1)
    DataIn(end+1:size(mask,1), :, :) = NaN; %add some NaNs
end
if dSize(2) < size(mask,2)
    DataIn(:, end+1:size(mask,2), :) = NaN; %add some NaNs
end
   
DataIn = reshape(DataIn,[numel(mask),prod(dSize(ndims(mask)+1:end))]); %merge x and y dimension based on mask size and remaining dimensions.
mask = mask(:); %reshape mask to vector
DataIn(mask,:) = [];
DataIn = reshape(DataIn,[size(DataIn,1),dSize(ndims(mask)+1:end)]);

% convert back to original size
dSize = size(DataIn); %size of input matrixmSize = size(mask);
mask = mask(:); %reshape mask to vector
DataOut = NaN([numel(mask) dSize(2:end)],dType); %pre-allocate new matrix
DataOut(~mask,:) = reshape(DataIn,sum(~mask),[]);
DataOut = reshape(DataOut,[mSize dSize(2:end)]);

end