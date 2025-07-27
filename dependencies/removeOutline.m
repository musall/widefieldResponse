function DataIn = removeOutline(DataIn,cropSize,cVal)
% short code to crop the outer edge of an image that was created using
% 'arrayShrink'. This means that there are NaNs in the background that can
% be used to create a mask. cVal can be used to add additional values to
% the mask.
% Usage: DataOut = removeOutline(DataIn,cropSize,cVal)
% -DataIn: Input image
% -cropSize: number of pixels to remove from the edge
% -cVal: background values behind the image. Default is NaNs

dSize = size(DataIn);
mask = isnan(DataIn(1:dSize(1),1:dSize(2)));
DataIn = reshape(DataIn,dSize(1),dSize(2),[]);

for x = 1 : size(DataIn,3)
    if exist('cVal','var') && ~isempty(cVal)
        cMask = mask | DataIn(:,:,x) == cVal;
    else
        cMask = mask;
    end
    cMask = imdilate(cMask,strel('square',cropSize));
    if length(dSize) == 3
        DataIn(:,:,x) = arrayCrop(DataIn(:,:,x),cMask);
    else
        DataIn = arrayCrop(DataIn,cMask);
    end
end
DataIn = reshape(DataIn,dSize);
