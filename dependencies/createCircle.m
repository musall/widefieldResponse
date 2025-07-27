function mask = createCircle(imSize, center, radius)
% Create a logical image of a circle with specified diameter, center, and 
% image size. This is useful to create a circular mask when analyzing
% imaging data. 
% Usage: createCircle(imSize, center, radius)
% - imSize: Size of the output image, as given by imSize = size(image)
% - center: X- and Y- coordinates for the circle center.
% - diamter: Diamter of the circle in pixels.

% First create the image.
[columnsInImage, rowsInImage] = meshgrid(1:imSize(2), 1:imSize(1));

% Next create the circle in the image.
mask = (rowsInImage - center(1)).^2 + (columnsInImage - center(2)).^2 <= radius.^2;