function [areaMask, areaLabels] = getAllenAreas(mask, minSize, splitM2)
%code to isolate individual areas from the allen atlas. Returns logical
%masks for each area on both hemisspheres and the correpsonding label.

load('allenDorsalMapSM.mat', 'dorsalMaps')
allenMask = dorsalMaps.allenMask;

if ~exist('mask','var') || isempty(mask)
    mask = false(size(allenMask)); %reject pixels within this mask
end

if ~exist('splitM2','var') || isempty(splitM2)
    splitM2 = false; %dont split M2 if this flag is false
end

% use allenMask if given as string
if strcmpi(mask, 'allenMask')
    mask = allenMask;
end

if ~exist('minSize','var') || isempty(minSize)
    minSize = 0; %minimum area size
end

areaMask = {}; areaLabels = {};
leftIdx = find(ismember(dorsalMaps.sidesSplit,'L'));
rightIdx = find(ismember(dorsalMaps.sidesSplit,'R'));
Cnt = 0;
for iAreas = 1 : length(leftIdx)
    cIdx = [leftIdx(iAreas) rightIdx(iAreas)];
    cOutline = dorsalMaps.edgeOutlineSplit{leftIdx(iAreas)};
    cMask = poly2mask(cOutline(:,2), cOutline(:,1),size(allenMask,1),size(allenMask,2));
    cOutline = dorsalMaps.edgeOutlineSplit{rightIdx(iAreas)};
    cMask = cMask | poly2mask(cOutline(:,2), cOutline(:,1),size(allenMask,1),size(allenMask,2));
    cMask = arrayCrop(cMask,mask);
    
    if strcmpi(dorsalMaps.labelsSplit{cIdx(1)}, 'MOs') && splitM2
        Cnt = Cnt + 1;
        circleMask = createCircle(size(cMask), size(cMask)/2, 100);
        areaMask{Cnt} = cMask == 1 & ~circleMask;
        areaLabels{Cnt} = 'M2a';
        
        Cnt = Cnt + 1;
        areaMask{Cnt} = cMask  == 1& circleMask;
        areaLabels{Cnt} = 'M2p';
    else
        if nansum(cMask(:)) > minSize
            Cnt = Cnt + 1;
            areaMask{Cnt} = cMask == 1;
            areaLabels{Cnt} = dorsalMaps.labelsSplit{cIdx(1)};
        end
    end
end