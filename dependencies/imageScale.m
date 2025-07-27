function [mapImg, cRange] = imageScale(cMap, cRange, twoSided, colors)
% quick command to plot image and set NaNs to be transparent. cRange
% defines a symmetric color range, based on top 'cRange' percentile in the image.
% 'twoSided' determines whether color range should be positive and negative
% (default) or only positive.
% usage: [mapImg, cRange] = imageScale(cMap, cRange, twoSided)

if ~exist('cRange', 'var') || isempty(cRange)
    cRange = abs(prctile(cMap(:),97.5));
end

if ~exist('twoSided', 'var') || isempty(twoSided)
    twoSided = true;
end

if ~exist('colors', 'var')
    colors = colormap_blueblackred;
end

if twoSided && isscalar(cRange)
    cRange = [-cRange cRange];
elseif isscalar(cRange)
    cRange = [prctile(cMap(:),2.5) prctile(cMap(:),97.5)];
end

mapImg = imshow(squeeze(cMap),cRange, 'InitialMagnification', 'fit');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
if ~strcmp(colors, "gray")
    colormap(mapImg.Parent, colors);
end




% nice colormap
function map = colormap_blueblackred(m)
if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

map=[1 1 0; 1 0.96 0; 1 0.92 0; 1 0.88 0; 1 0.84 0; 1 0.80 0; 1 0.76 0; 1 0.72 0; 1 0.68 0; 1 0.64 0; 1 0.60 0; 1 0.56 0; 1 0.52 0; 1 0.48 0; 1 0.44 0; 1 0.40 0;  ...
    1 0.36 0; 1 0.32 0; 1 0.28 0; 1 0.24 0; 1 0.20 0; 1 0.16 0; 1 0.12 0; 1 0.08 0; 1 0.04 0;
    1 0 0; 0.96 0 0; 0.92 0 0; 0.88 0 0; 0.84 0 0; 0.80 0 0; 0.76 0 0; 0.72 0 0; 0.68 0 0; 0.64 0 0; 0.60 0 0; 0.56 0 0; 0.52 0 0; 0.48 0 0; 0.44 0 0; 0.40 0 0;  ...
    0.36 0 0; 0.32 0 0; 0.28 0 0; 0.24 0 0; 0.20 0 0; 0.16 0 0; 0.12 0 0; 0.08 0 0; 0.04 0 0; 0 0 0;                                   ...
    0 0 0.04;  0 0 0.08; 0 0 0.12; 0 0 0.16; 0 0 0.20; 0 0 0.24; 0 0 0.28; 0 0 0.32; 0 0 0.36; 0 0 0.40; 0 0 0.44; 0 0 0.48; 0 0 0.52; ...
    0 0 0.56; 0 0 0.60; 0 0 0.64; 0 0 0.68; 0 0 0.72; 0 0 0.76; 0 0 0.80; 0 0 0.84; 0 0 0.88; 0 0 0.92; 0 0 0.96; 0 0 1; ...
    0 0.04 1;  0 0.08 1; 0 0.12 1; 0 0.16 1; 0 0.20 1; 0 0.24 1; 0 0.28 1; 0 0.32 1; 0 0.36 1; 0 0.40 1; 0 0.44 1; 0 0.48 1; 0 0.52 1; ...
    0 0.56 1; 0 0.60 1; 0 0.64 1; 0 0.68 1; 0 0.72 1; 0 0.76 1; 0 0.80 1; 0 0.84 1; 0 0.88 1; 0 0.92 1; 0 0.96 1; 0 1 1];

map = imresize(map,[m,3]);
map = map - min(map(:));
map = map ./ max(map(:));
map = map(end:-1:1,:);