function oLines = cc_plotAllenOutline(ax,side,allRegions, lineWidths)
%function to add area outlines from allen CCF to current image

load('allenDorsalMapSM.mat')
if ~exist('ax','var') || isempty(ax)
    ax = gca;
end

cIdx = true(1,length(dorsalMaps.edgeOutlineSplitRed));
if exist('side','var')
    if strcmpi(side,'L') || strcmpi(side,'R')
        cIdx = ismember(dorsalMaps.sidesSplitRed,side)';
    end
end

if ~exist('allRegions','var')
    allRegions = false;
end

if ~exist('lineWidths','var')
    lineWidths = 0.1;
end

if ishold(ax)
    checker = true;
else
    hold(ax,'on'); checker = false;
end

%check image size
shiftMidline = 0;
for x = 1 : length(ax.Children)
    if contains(class(ax.Children(x)),'Image')
        
        allenMask = dorsalMaps.allenMask; %allenMask
        cImg = ax.Children(x).CData; %current image
        imgSize = size(cImg); %size of current image
        lineScale = min((size(allenMask) ./ imgSize(1:2))); %scale lines to match size of the image
        
        % check if there is an offset in the midline if camera image is
        % square instead of rectangular (as expected by allenMask)
        if size(cImg,1) == size(cImg,2)
            cMidline = (size(allenMask, 2) / 2); %current midline
            tMidline = (size(cImg, 2) / 2); %target midline
            shiftMidline = cMidline - tMidline; %how much does the mask need to be shifted?
        else
            shiftMidline = 0;
        end

    end
end

for x = find(cIdx)
    if allRegions
        oLines(x) = plot(ax, dorsalMaps.edgeOutlineSplit{x}(:,2)./lineScale - shiftMidline, dorsalMaps.edgeOutlineSplit{x}(:,1)./lineScale, 'w', 'LineWidth', lineWidths);
    else
        oLines(x) = plot(ax, dorsalMaps.edgeOutlineSplitRed{x}(:,2)./lineScale - shiftMidline, dorsalMaps.edgeOutlineSplitRed{x}(:,1)./lineScale, 'w', 'LineWidth', lineWidths);
    end
end

if ~checker
    hold(ax,'off');
end