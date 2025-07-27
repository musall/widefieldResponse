%% data path
recs = {
    'F:\widefieldResponse_data\2897\Rec1'
    'F:\widefieldResponse_data\2897\Rec2'
    'F:\widefieldResponse_data\2897\Rec3'
    'F:\widefieldResponse_data\2908\Rec1'
    'F:\widefieldResponse_data\2908\Rec2'
    'F:\widefieldResponse_data\2908\Rec3'
    'F:\widefieldResponse_data\2909\Rec1'
    'F:\widefieldResponse_data\2909\Rec2'
    'F:\widefieldResponse_data\2909\Rec3'
    'F:\widefieldResponse_data\2910\Rec1'
    'F:\widefieldResponse_data\2910\Rec2'
    'F:\widefieldResponse_data\2910\Rec3'
    };

%% settings
opts2.firstDur = 10;
opts2.secondDur = 60;
opts2.postStimDur = 300;
opts2.baselineDur = 240;
opts2.sRate = 15.1515;
opts2.areaCoords = [335,180; 370,260; 375,415];
opts2.areaSize = 40; %size of ROIs for traces
opts2.areaNames = {'Frontal', 'Motor', 'Parietal'};
opts2.sharePath = 'F:\wr_results\';
opts2.crossLags = 10; %amount of shift for cross correlation in seconds
opts2.crossCorrTimeRange = 300; %time range for data in the crosscorrelation analysis (pre and post stim)

% for recovery analysis
recoveryWindow = 30; %window size used for recovery and amplitude analysis
recoveryWindowShift = 5; %shift of window used for recovery analysis
troughSamples = 150;
binWindow = 10; %winSize for binned maps
troughRange = 90; %time range to check for trough before and after stim

% other things
winRange = [-opts2.baselineDur opts2.postStimDur]; %range of cutout in seconds
gapWidth = round(5 * opts2.sRate); %nr of datapoints for nan gap during stimulation
earlyColor = [0.25 0.25 0.25];
lateColor = [0 0.5 0.1];
wfMoveWin = -15:30; %range of frames around movement events

load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask; %allenMask

%% load widefield / motion data from different recordings
nrRecs = size(recs,1);
baseFrames = ceil(opts2.baselineDur * opts2.sRate);
earlyStimFrames = ceil((opts2.firstDur) * opts2.sRate);
lateStimFrames = ceil((opts2.secondDur) * opts2.sRate);
stimFrames = ceil(opts2.postStimDur * opts2.sRate);

overviewFig = figure;
t = tiledlayout(4,3);
t.TileSpacing = 'compact'; % 'compact' or 'none'
t.Padding = 'compact'; % 'compact' or 'none'
allU = cell(1,nrRecs);
allV = cell(1,nrRecs);
whiskerME = cell(1,nrRecs);
stimMotion = nan(3, nrRecs);
for iRecs = 1 : nrRecs %current recording
    
    % get some information for current recording
    basePath = recs{iRecs};
    [~, recName] = (fileparts(basePath));
    [~, cAnimal] = fileparts(fileparts(basePath));
    
    %% load vessel image and align to Allen
    clear transParams opts
    load(fullfile(basePath, 'blueAvg.mat'), 'blueAvg');
    load(fullfile(basePath, 'opts2.mat'), 'opts');
    transParams = opts.transParams;
    vesselPic = alignAllenTransIm(nanmean(blueAvg,3), transParams);
    
    nexttile
    imshow(mat2gray(vesselPic))
    cc_plotAllenOutline;
    title([cAnimal '; ' recName]);
    drawnow;
    
    %% get functional data
    imgFile = fullfile(basePath, 'nVc.mat');
    load(imgFile, 'nVc', 'U', 'preIdx', 'postIdx', 'preStimDur');
    
    % align to Allen
    U = alignAllenTransIm(U, transParams);
    U = rateDisc_removeOutline(U, 10);
    U(isnan(U)) = 0;
    allU{iRecs} = arrayCrop(U, allenMask);
    
    % get stimulus response
    baselineIdx = preStimDur - baseFrames + 1 : preStimDur;
    stimIdx = preStimDur + (1 : stimFrames);
    allV{iRecs} = nVc(:,[baselineIdx, stimIdx]);
    
    %% get video data
    cFile = dir([basePath, filesep, '*cam1*frameTimes*']);
    frameTimes = fullfile(cFile.folder, cFile.name);
    load(frameTimes, 'frameInfo');
    frameTimes = frameInfo(:,2) - frameInfo(1,2);
    
    load([basePath, filesep, 'motionEnergy.mat'], 'vidTrigs', 'vidTimes');
    load([basePath, filesep, 'faceVars.mat'], 'whiskerMotion');
    whiskerMotion = mV_vidResamp(whiskerMotion, frameTimes', 30)';
    
    cTrig = find(vidTimes{1} > vidTrigs{1}(1), 1) ;
    baselineIdx = cTrig - ceil(30*opts2.baselineDur) : cTrig;
    stimIdx = cTrig : cTrig + ceil(30*opts2.postStimDur);
    
    % compute motion energy before, during, and after stimulus
    whiskerME{iRecs} = whiskerMotion([baselineIdx, stimIdx]);
    
    loopReporter(iRecs, nrRecs, 10)
end

%% show activity maps for equal bins afters
timeAxis = (1 : stimFrames + baseFrames + gapWidth) / opts2.sRate;
timeAxis = timeAxis - baseFrames / opts2.sRate;
areaSize = opts2.areaSize;
cMask = allenMask;
cMask(:, 1: size(allenMask, 2) / 2) = true; %only right HS

binRange = -(3*binWindow) : binWindow : opts2.postStimDur;
nrBins = length(binRange) - 1;
binData = nan(sum(~cMask(:)), nrBins, nrRecs, 'single');
for iRecs = 1 : nrRecs
    cU = arrayShrink(allU{iRecs}, cMask, 'merge');
    baseline = cU * nanmean(allV{iRecs}(:, 1 : baseFrames),2);
    for iBins = 1 : nrBins
        cIdx = timeAxis > binRange(iBins) & timeAxis <= binRange(iBins+1);
        cV = nanmean(allV{iRecs}(:, cIdx),2) - median(allV{iRecs}(:, 1 : baseFrames),2);
        binData(:, iBins, iRecs) = cU * cV;
    end
end
binData = arrayShrink(binData, cMask, 'split');

%% make figure - average bins
h = figure('name', 'Average response maps');
t = tiledlayout(3, ceil(nrBins/3));
t.TileSpacing = 'compact'; % 'compact' or 'none'
t.Padding = 'compact'; % 'compact' or 'none'
cRange = 0.025;

for iBins = 1 : nrBins
    nexttile;
    imageScale(nanmean(binData(:,:,iBins,:), 4), cRange, true); hold on;
    for x = 1 : size(opts2.areaCoords, 1)
        rectangle('Position',[opts2.areaCoords(x,:)-areaSize/2, repmat(areaSize,1,2)],'linewidth',2, 'Curvature',1, 'edgecolor', 'w', 'parent', gca);
    end
    cc_plotAllenOutline(gca, 'R');
    title(sprintf('%i - %i seconds', binRange(iBins), binRange(iBins+1)));
    xlim([round(size(allenMask, 2) / 2), size(allenMask, 2)]);
    drawnow;
end

%% single rec examples
iRec = 2;
[~, recName] = (fileparts(recs{iRecs}));
[~, cAnimal] = fileparts(fileparts(recs{iRecs}));

h = figure('name',[cAnimal '-' recName], 'renderer', 'painters');
t = tiledlayout(3, ceil(nrBins/3));
t.TileSpacing = 'compact'; % 'compact' or 'none'
t.Padding = 'compact'; % 'compact' or 'none'
cRange = 0.04;

for iBins = 1 : nrBins
    nexttile;
    imageScale(binData(:,:,iBins,iRec), cRange, true); hold on;
    for x = 1 : size(opts2.areaCoords, 1)
        rectangle('Position',[opts2.areaCoords(x,:)-areaSize/2, repmat(areaSize,1,2)],'linewidth',2, 'Curvature',1, 'edgecolor', 'w', 'parent', gca);
    end
    cc_plotAllenOutline(gca, 'R');
    title(sprintf('%i - %i seconds', binRange(iBins), binRange(iBins+1)));
    xlim([round(size(allenMask, 2) / 2), size(allenMask, 2)]);
    drawnow;
end

%% show activity maps for early and late responses
cRange = 0.02;
areaSize = opts2.areaSize;
cData = cell(1,2);
for iRecs = 1 : nrRecs
    baseline = svdFrameReconstruct(allU{iRecs}, nanmean(allV{iRecs}(:, 1 : baseFrames),2));
    
    %early response
    stim1 = svdFrameReconstruct(allU{iRecs}, nanmean(allV{iRecs}(:, baseFrames : baseFrames + earlyStimFrames),2)) - baseline;
    cData{1} = runMean(cData{1}, stim1, iRecs);
    
    % late response
    stim2 = svdFrameReconstruct(allU{iRecs}, nanmean(allV{iRecs}(:, baseFrames + (earlyStimFrames : lateStimFrames)),2)) - baseline;
    cData{2} = runMean(cData{2}, stim2, iRecs);
end

cTitles = {['Early response: First ' int2str(opts2.firstDur) ' seconds'], ...
    ['Late response: Later ' int2str(opts2.secondDur) ' seconds']};
cColors = {'k', 'w'};

figure
for iStims = 1 : 2
    subplot(1,2,iStims)
    cImg = imageScale(cData{iStims});
    colormap(cImg.Parent, colormap_blueblackred);
    caxis([-cRange cRange]);
    cc_plotAllenOutline(gca);
    title(cTitles{iStims});
    hold on;
    for x = 1 : size(opts2.areaCoords, 1)
        rectangle('Position',[opts2.areaCoords(x,:)-areaSize/2, repmat(areaSize,1,2)],'linewidth',2, 'Curvature',1, 'edgecolor', cColors{iStims}, 'parent',cImg.Parent);
    end
end

%% get activity for each araea
[areaMask, areaLabels] = getAllenAreas('allenMask', 10, true);
nrAreas = size(opts2.areaCoords,1);
dSize = size(allU{1});

h = figure('renderer', 'painters');
h.Position = [681    42   560   954];
t = tiledlayout(nrAreas,1);
t.TileSpacing = 'compact'; % 'compact' or 'none'
t.Padding = 'compact'; % 'compact' or 'none'
areaAmp = nan(nrAreas, nrRecs, 2, 'single');
traceData = cell(1, nrAreas);
for iAreas = 1 : nrAreas
    
    usePixels = createCircle(dSize(1:2), fliplr(opts2.areaCoords(iAreas,:)), opts2.areaSize);
    traceData{iAreas} = nan(stimFrames + baseFrames, nrRecs, 'single');
    for iRecs = 1 : nrRecs
        temp = reshape(allU{iRecs}, [], dSize(3));
        cTrace = nanmean(temp(usePixels,:) * allV{iRecs},1);
        traceData{iAreas}(:, iRecs) = cTrace - prctile(cTrace(1:baseFrames), 50);
    end
    
    %some quantification
    areaAmp(iAreas,:,1) = nanmean(traceData{iAreas}(baseFrames : baseFrames + earlyStimFrames, :));
    areaAmp(iAreas,:,2) = nanmean(traceData{iAreas}(baseFrames + earlyStimFrames : baseFrames + stimFrames, :));
    
    % show traces
    traceData{iAreas} = insertColumns(traceData{iAreas}', baseFrames+1, nan, gapWidth);
    
    nexttile
    arrayPlot(timeAxis, traceData{iAreas}', 'k')
    ylim([-0.1 0.15]);
    nhline(0, '--k');
    ylabel('dF/F');
    xlabel('time after stimulus onset (s)');
    title(['Area ' num2str(iAreas)]);
    nvline([0, 5], '--k');
    xlim([-opts2.baselineDur, opts2.postStimDur]);
    ax = gca;
    ax.FontSize = 14;
    ax.TickDir = 'out';
    set(gca,'box','off')
    
end

%% compute cross-corelation between areas before and after stim
% loop over areas and check if their correlations are time-shifted or if peak correlations occur without delay
shiftLength = round(opts2.sRate*opts2.crossLags);
corrPre = nan(nrAreas, nrAreas, nrRecs, 'single');
corrPost = nan(nrAreas, nrAreas, nrRecs, 'single');
crossCorrSpectraPre = nan(nrAreas, nrAreas, shiftLength*2+1, nrRecs, 'single');
crossCorrSpectraPost = nan(nrAreas, nrAreas, shiftLength*2+1, nrRecs, 'single');
crossCorrMaxPre = nan(nrAreas, nrAreas, nrRecs, 'single');
crossCorrMaxPost = nan(nrAreas, nrAreas, nrRecs, 'single');
crossCorrLagPre = nan(nrAreas, nrAreas, nrRecs, 'single');
crossCorrLagPost = nan(nrAreas, nrAreas, nrRecs, 'single');
for iRecs = 1 : nrRecs
    for x = 1 : nrAreas
        
        crossCorrFrames = round(opts2.crossCorrTimeRange .* opts2.sRate) + 1; %number of frames to be used in each period
        preStimStart = max([1, baseFrames - crossCorrFrames]);
        postStimEnd = min([size(traceData{1},2), baseFrames + crossCorrFrames + gapWidth]);
        
        preData = traceData{x}(iRecs, preStimStart:baseFrames)';
        postData = traceData{x}(iRecs, baseFrames + gapWidth + 1 : postStimEnd)';
        
        for y = 1 : nrAreas
            
            %preStim cross correlation
            cData = traceData{y}(iRecs, preStimStart:baseFrames)';
            cSpectra = crosscorr(preData, cData, 'NumLags', shiftLength);
            corrPre(x,y, iRecs) = cSpectra(shiftLength+1); %correlation at 0 shift. This is the same as regular correlation.
            [crossCorrMaxPre(x,y, iRecs), crossCorrLagPre(x,y, iRecs)] = max(cSpectra);
            crossCorrSpectraPre(x,y,:, iRecs) = cSpectra;
            
            %postStim cross correlation
            cData = traceData{y}(iRecs, baseFrames + gapWidth + 1 : postStimEnd)';
            cSpectra = crosscorr(postData, cData, 'NumLags', shiftLength);
            corrPost(x,y, iRecs) = cSpectra(shiftLength+1); %correlation at 0 shift. This is the same as regular correlation.
            [crossCorrMaxPost(x,y, iRecs), crossCorrLagPost(x,y, iRecs)] = max(cSpectra);
            crossCorrSpectraPost(x,y,:, iRecs) = cSpectra;
            
        end
    end
end

%subtract the spectra center to compute symmetric relative time lag
crossCorrLagPre = crossCorrLagPre - shiftLength - 1;
crossCorrLagPost = crossCorrLagPost - shiftLength - 1;

%% show correlations results in a figure
h = figure('renderer', 'painters');
subplot(2,3,[1 4]);
xRange = (-152 : shiftLength)./opts2.sRate;
arrayPlot(xRange, squeeze(crossCorrSpectraPost(1,3,:,:)), 'k'); axis square
xlim([-shiftLength shiftLength]./opts2.sRate)
grid on;
xlabel('Time shift (s)')
ylabel('Correlation coefficient')
title('Cross-correlation: Frontal-PPC')
ax = gca;
ax.YTick = -0.5:0.25:1;
ax.XTick = -10:5:10;
niceFigure;

cTitles = {'', 'Max corr - PreStim', 'Corr - PreStim', '', 'Max corr - PostStim', 'Corr - PostStim'};
for x = [2, 3, 5, 6]
    subplot(2,3,x)
    if x == 2
        imagesc(nanmean(crossCorrMaxPre,3));
    elseif x == 3
        imagesc(nanmean(corrPre,3));
    elseif x == 5
        imagesc(nanmean(crossCorrMaxPost,3));
    elseif x == 6
        imagesc(nanmean(corrPost,3));
    end
    cb = colorbar; axis square
    caxis([0.7 1]);
    ylabel(cb, 'Correlation coefficient', 'FontSize', 14, 'Rotation', 270);
    cb.Label.Position(1) = 4;
    ax = gca;
    colormap(ax, hot(256));
    title(cTitles{x});
    ax.XTick = 1:nrAreas;
    ax.YTick = 1:nrAreas;
    ax.XTickLabel = opts2.areaNames;
    ax.YTickLabel = opts2.areaNames;
    ax.TickLength = [0 0];
end

%% show response amplitudes
figure
t = tiledlayout(nrAreas,1);
t.TileSpacing = 'compact'; % 'compact' or 'none'
t.Padding = 'compact'; % 'compact' or 'none'

clear h
for iAreas = 1 : nrAreas
    nexttile
    h(1) = Violin(areaAmp(iAreas,:,1).*100, 0, 'ViolinColor', earlyColor, 'EdgeColor', earlyColor, 'BoxColor', [0 0 0], 'Bandwidth', 0.75);
    h(2) = Violin(areaAmp(iAreas,:,2).*100, 1, 'ViolinColor', lateColor, 'EdgeColor', lateColor, 'BoxColor', [0 0 0], 'Bandwidth', 0.75);
    ylim([-6 6])
    xlim([-1 2])
    axis square
    ax = h(1).ViolinPlot.Parent;
    ax.XTick = 0:1;
    ax.XTickLabel = {'Early', 'Late'};
    ylabel('dF/F (%)');
    title(['Area ' num2str(iAreas)]);
    nhline(0, 'k--')
    niceFigure;
end

%% example traces individual recording
plotRange = [-30 90];
preStimAmp = nan(nrAreas, nrRecs);
preStimMedian = nan(nrAreas, nrRecs);
latePostAmp = nan(nrAreas, nrRecs);
earlyPreAmp = nan(nrAreas, nrRecs);
earlyPostAmp = nan(nrAreas, nrRecs);
preStimTime = nan(nrAreas, nrRecs);
postStimTime = nan(nrAreas, nrRecs);
recoveryTime = nan(nrAreas, nrRecs);
recoveryAmp = nan(nrAreas, nrRecs);
whiskVid = cell(1, nrRecs);
for iRecs = 1 : nrRecs
    
    % get some information for current recording
    basePath = recs{iRecs};
    [~, recName] = (fileparts(basePath));
    [~, cAnimal] = fileparts(fileparts(basePath));
    targChan = 1;
    
    clear eData
    h5Files = dir(fullfile(basePath, '*.h5'));
    [ephysData, analogData, digitalData, sRate] = load_MCSh5file_SM(basePath,h5Files(1).name);
    stimDur = digitalData.events{1}{3}(1) - digitalData.events{1}{2}(1); %stimulus duration in seconds
    stimOn = digitalData.events{1}{2}(1); %trigger onset in seconds
    useIdx = (round(stimOn * sRate) + (winRange(1) * sRate : winRange(2) * sRate -1))';
    useIdx(useIdx > size(ephysData.samples,2)) = [];
    eData = ephysData.samples(targChan,useIdx)';
    eData = eData - mean(eData(round((abs(winRange(1))-10) * sRate):round(abs(winRange(1)) * sRate)));
    xRange = winRange(1) : 1/sRate : winRange(2)-1/sRate;
    
    if iRecs == 1 %show example traces
        h = figure('renderer', 'painters');
        h.Position = [681    42   560   954];
        t = tiledlayout(nrAreas+2,1);
        t.TileSpacing = 'compact'; % 'compact' or 'none'
        t.Padding = 'compact'; % 'compact' or 'none'
        
        % ephys data
        nexttile
        plot(xRange(1:length(eData)), eData, 'k');
        xlabel('time from trigger (s)')
        ylabel('voltage')
        xlim(winRange);
        nhline(0, '--k');
        nvline([0 5], '--k');
        title(basePath);
        ax = gca;
        ax.XTickLabel = [];
        ax.TickDir = 'out';
        ax.XColor = 'w';
        
        % motion energy
        nexttile
        xRangeME = winRange(1) : 1/30 : winRange(2)-1/30;
        cData = smooth(whiskerME{iRecs}(1:length(xRangeME)),30);
        plot(xRangeME, cData, 'k');
        xlabel('time from trigger (s)');
        ylabel('Motion energy [AU]');
        ylim([0 max(cData)+1]);
        xlim(winRange);
        nhline(0, '--k');
        nvline([0 5], '--k');
        title('Motion energy');
        ax = gca;
        ax.XTickLabel = [];
        ax.TickDir = 'out';
        ax.XColor = 'w';
    end
    
    %% check for video ME events
    temp = smoothCol(whiskerME{iRecs}(1:length(xRangeME)),2, 10, 'box');
    temp = temp - prctile(temp, 50);
    temp = temp ./ std(temp);
    eTrigs = xRangeME(diff(temp > 2.5) == 1);
    
    % find events in the widefield
    wfTrigs = nan(1, length(eTrigs));
    for iTrigs = 1 : length(eTrigs)
        [~, wfTrigs(iTrigs)] = min(abs(timeAxis - eTrigs(iTrigs)));
    end
    wfTrigs = wfTrigs(wfTrigs > abs(wfMoveWin(1)));
    wfTrigs = wfTrigs(wfTrigs + wfMoveWin(end) < length(timeAxis));
    baseTrigs = wfTrigs(timeAxis(wfTrigs) < 0)  + wfMoveWin';
    postTrigs = (wfTrigs(timeAxis(wfTrigs) > 0)-18)  + wfMoveWin';
    
    cV = allV{iRecs};
    cV = insertColumns(cV, baseFrames+1, nan, gapWidth);
    baseTrigV = cV(:, baseTrigs(:));
    baseTrigV = nanmedian(reshape(baseTrigV, [size(baseTrigV,1), size(baseTrigs)]),3);
    baseTrigV = baseTrigV - baseTrigV(:,1);
    baseTrigV = svdFrameReconstruct(allU{iRecs}, baseTrigV);
    
    postTrigV = cV(:, postTrigs(:));
    postTrigV = nanmedian(reshape(postTrigV, [size(postTrigV,1), size(postTrigs)]),3);
    postTrigV = postTrigV - postTrigV(:,1);
    postTrigV = svdFrameReconstruct(allU{iRecs}, postTrigV);
    whiskVid{iRecs} = cat(4, baseTrigV, postTrigV);
    
    for iAreas = 1 : nrAreas
        
        % get dF/F traces
        cData = traceData{iAreas}(iRecs,:);
        smoothData = cData;
        smoothData(1:find(isnan(cData),1)-1) = smoothCol(cData(1:find(isnan(cData),1)-1),2,451);
        smoothData(find(isnan(cData),1, 'last')+1:end) = smoothCol(cData(find(isnan(cData),1, 'last')+1:end),2,451);
        
        cData(1:find(isnan(cData),1)-1) = smoothCol(cData(1:find(isnan(cData),1)-1),2,15);
        cData(find(isnan(cData),1, 'last')+1:end) = smoothCol(cData(find(isnan(cData),1, 'last')+1:end),2,15);
        
        if iRecs == 1
            nexttile
            plot(timeAxis, cData, 'k');  hold on;
            plot(timeAxis, smoothData, 'color', [0.5 0.5 0.5], 'linewidth', 2);
            
            ylabel('dF/F');
            xlabel('time after stimulus onset (s)');
            title(['Area ' num2str(iAreas)]);
            xlim(winRange);
            nhline(0, '--k');
            nvline([0 5], '--k');
            ax = gca;
            ax.FontSize = 14;
            ax.TickDir = 'out';
            set(gca,'box','off')
            ax = gca;
            if iAreas ~= nrAreas
                ax.XColor = 'w';
            end
        end
        
        % get lowest values in baseline
        preStimIdx = timeAxis < 0 & timeAxis > -troughRange;
        preStimTimeAxis = timeAxis(preStimIdx);
        preStimData = smoothData(preStimIdx); %smoothed poststim data
        lowPreStimIdx = find(preStimData < prctile(preStimData, 5));
        preStimTime(iAreas, iRecs) = median(preStimTimeAxis(lowPreStimIdx)); %median time of 5% smalles values
        cTimeIdx = preStimTimeAxis > preStimTime(iAreas, iRecs) - (recoveryWindow/2) & preStimTimeAxis < preStimTime(iAreas, iRecs) + (recoveryWindow/2);
        preStimAmp(iAreas, iRecs) = nanmedian(preStimData(cTimeIdx))*100; %median amplitude of 5% smalles values as percentage
        preStimMedian(iAreas, iRecs) = nanmedian(preStimData)*100; %median of baseline as percentage
        
        % get median for early response
        earlyStimIdx = timeAxis > (gapWidth / opts2.sRate) & timeAxis < opts2.firstDur + gapWidth / opts2.sRate;
        earlyPostAmp(iAreas, iRecs) = nanmedian(cData(earlyStimIdx)*100); %early poststim data
        
        earlyStimIdx = timeAxis < 0 & timeAxis > -opts2.firstDur;
        earlyPreAmp(iAreas, iRecs) = nanmedian(cData(earlyStimIdx)*100); %early prestim data
        
        % same thing for poststim but also find time points of lowest
        % values to identify time of trough
        postStimIdx = timeAxis > 0 & timeAxis < troughRange;
        postStimTimeAxis = timeAxis(postStimIdx);
        postStimData = smoothData(postStimIdx); %smoothed poststim data
        lowPostStimIdx = find(postStimData < prctile(postStimData, 5));
        postStimTime(iAreas, iRecs) = median(postStimTimeAxis(lowPostStimIdx)); %median time of 5% smalles values
        cTimeIdx = postStimTimeAxis > postStimTime(iAreas, iRecs) - (recoveryWindow/2) & postStimTimeAxis < postStimTime(iAreas, iRecs) + (recoveryWindow/2);
        latePostAmp(iAreas, iRecs) = nanmedian(postStimData(cTimeIdx))*100; %median amplitude of 5% smalles values as percentage
        
        % compute recovery time after the through
        winSteps = postStimTime(iAreas, iRecs)  : recoveryWindowShift : postStimTimeAxis(end)-recoveryWindow;
        cAmpData = nan(1, length(winSteps));
        Cnt = 0;
        postStimData = smoothData(postStimIdx); %smoothed poststim data
        for iWindows = winSteps
            
            Cnt = Cnt + 1;
            cTimeIdx = postStimTimeAxis > iWindows - (recoveryWindow/2) & postStimTimeAxis < iWindows + (recoveryWindow/2);
            cWinData = postStimData(cTimeIdx); %data for current window
            lowIdx = cWinData < prctile(cWinData, 5); %5% of smallest values in the current window
            cAmpData(Cnt) = nanmedian(cWinData)*100; %median of 5% smalles values as percentage in current window
            
            if cAmpData(Cnt) > preStimAmp(iAreas, iRecs) %recovery point is reached when min values in current window are above min values in the baseline
                recoveryTime(iAreas, iRecs) = iWindows;
                recoveryAmp(iAreas, iRecs) = cAmpData(Cnt);
                break
            end
        end
        
        % find trough before the stim to compare amplitude
        winSteps = postStimTime(iAreas, iRecs)  : recoveryWindowShift : winRange(end)-recoveryWindow;
        cAmpData = nan(1, length(winSteps));
        Cnt = 0;
        postStimData = smoothData(postStimIdx); %smoothed poststim data
        for iWindows = winSteps
            
            Cnt = Cnt + 1;
            cTimeIdx = postStimTimeAxis > iWindows - (recoveryWindow/2) & postStimTimeAxis < iWindows + (recoveryWindow/2);
            cWinData = postStimData(cTimeIdx); %data for current window
            lowIdx = cWinData < prctile(cWinData, 5); %5% of smallest values in the current window
            cAmpData(Cnt) = nanmedian(cWinData)*100; %median of 5% smalles values as percentage in current window
            
            if cAmpData(Cnt) > preStimAmp(iAreas, iRecs) %recovery point is reached when min values in current window are above min values in the baseline
                recoveryTime(iAreas, iRecs) = iWindows;
                recoveryAmp(iAreas, iRecs) = cAmpData(Cnt);
                break
            end
        end
        
        if iRecs == 1
            %add a vertical line to indicate the trough and recovery time
            nvline(preStimTime(iAreas, iRecs), 'g--');
            nvline(postStimTime(iAreas, iRecs), 'r--');
            nvline(recoveryTime(iAreas, iRecs), 'b--');
            drawnow;
        end
    end
end

%% motion triggered widefield result (needs data from all individual recordings in the loop above)
% show images
useIdx = 25:40;
cData = cat(5, whiskVid{:}).*100;
moveFig = figure;
groupNames = {'Before stim', 'After stim'};
for x = 1 : 2
    subplot(1,3,x);
    temp = squeeze(nanmean(nanmean(cData(:,:,useIdx,x,:),3),5));
    cImg = imageScale(temp);
    title(groupNames{x});
    colorbar;
    cc_plotAllenOutline;
    caxis([-4 4])
end

% show traces
usePixels = createCircle(dSize(1:2), fliplr([395,315]), opts2.areaSize);

datSize = size(cData);
cData = reshape(cData,size(cData,1)*size(cData,2), []);
cTrace = nanmean(cData(usePixels,:),1);
cTrace = reshape(cTrace, datSize(3), datSize(4), datSize(5));

subplot(1,3,3); clear cLine
moveTimeAx = wfMoveWin/opts.sRate;
cLine(1) = stdshade(squeeze(cTrace(:,1,:))', 0.2, 'r', moveTimeAx);
cLine(2) = stdshade(squeeze(cTrace(:,2,:))', 0.2, 'b', moveTimeAx);
legend(cLine, groupNames, 'Location', 'northwest');
axis square;
niceFigure
grid on;
xlabel('time after move (s)');
ylabel('dF / F (%)');
nhline(0, 'k');

%% do some statistics for movement triggered WF response
[maxVal,maxTime] = max(cTrace, [], 1);
maxVal = squeeze(maxVal);
maxTime = squeeze(maxTime) ./ 30;

cFig = figure;
subplot(1,2,1);
nConds = 2;
for x = 1 : nConds
    Violin(maxVal(x,:), x, 'ViolinColor', [0 0 0], 'EdgeColor', [0 0 0], 'BoxColor', [0 0 0]);
end
ax = gca;
ax.XTick = 1:nConds;
ax.XTickLabel = {'Before', 'During'};
p = signrank(maxVal(1,:), maxVal(2,:));
title(['Max. motion response - p signrank test = ' num2str(p)])
niceFigure;
ylim([0 12]);
axis square;
ylabel('peak motion response - dF / F (%)');

subplot(1,2,2);
for x = 1 : nConds
    Violin(maxTime(x,:), x, 'ViolinColor', [0 0 0], 'EdgeColor', [0 0 0], 'BoxColor', [0 0 0]);
end
ax = gca;
ax.XTick = 1:nConds;
ax.XTickLabel = {'Before', 'During'};
p = signrank(maxTime(1,:), maxTime(2,:));
title(['Delay to max response - p signrank test = ' num2str(p)])
niceFigure;
ylim([0.8 1.6]);
axis square;
ylabel('Delay to peak response (s)');


%% show some statistical results
h = figure;
t = tiledlayout(1,4);
t.TileSpacing = 'compact'; % 'compact' or 'none'
t.Padding = 'compact'; % 'compact' or 'none'

% result for early response
clear h
nexttile
earlyAmpData = earlyPostAmp - earlyPreAmp;
xlim([0 4]);
nhline(0, '--k');
for iAreas = 1 : nrAreas
    cData = earlyAmpData(iAreas,:);
    h(iAreas) = Violin(cData, iAreas, 'ViolinColor', earlyColor, 'EdgeColor', earlyColor, 'BoxColor', [0 0 0], 'BandWidth', 0.5);
    pEarly(iAreas) = signrank(cData(~isnan(cData(:)))); %check against 0
    
    % give feedback
    disp('==================')
    fprintf('Early response for Area %i:\n', iAreas)
    fprintf('PostStim = %.2f %c %.2f percent \n'  , mean(cData), char(177), sem(cData));
    fprintf('pVal ranksum test: %f\n', signrank(cData));
    disp('==================')
end

ax = h(1).ViolinPlot.Parent;
ylim([round(min(earlyAmpData(:)))-0.2 round(max(earlyAmpData(:)))+1])
axis square
ax.XTick = [1 : nrAreas];
xlabel('Area IDs')
ylabel('amplitude (dF/F)');
title('Early response');
subtitle(sprintf('Amps: %.4f, %.4f, %.4f\n pVals = %.6f, %.6f, %.6f', median(earlyAmpData, 2)', pEarly));

% late response
troughAmp = latePostAmp - preStimAmp;
nexttile
xlim([0 4]);
nhline(0, '--k'); clear h
for iAreas = 1 : nrAreas
    cData = troughAmp(iAreas,:);
    h(iAreas) = Violin(cData, iAreas, 'ViolinColor', earlyColor, 'EdgeColor', earlyColor, 'BoxColor', [0 0 0]);
    pLate(iAreas) = signrank(cData); %check against 0
    areaStats{1}(iAreas,:) = cData; %keep data to save
    
    % give feedback
    disp('==================')
    fprintf('Late response for Area %i:\n', iAreas)
    fprintf('PostStim = %.2f %c %.2f percent \n'  , mean(cData), char(177), sem(cData));
    fprintf('pVal ranksum test: %f\n', signrank(cData));
    disp('==================')
end
ax = h(1).ViolinPlot.Parent;
ylim([round(min(troughAmp(:)))-0.2 round(max(troughAmp(:)))+1])
axis square
ax.XTick = [1 : nrAreas];
xlabel('Area IDs')
ylabel('amplitude (dF/F)');
title('Late response (trough amp)');
subtitle(sprintf('Amps: %.4f, %.4f, %.4f\n pVals = %.6f, %.6f, %.6f', median(troughAmp, 2)', pLate));

% trough time
nexttile
useIdx = max(postStimTime, [], 1) < 90 & ~any(isnan(recoveryTime), 1); %only consider recovery for recordings with at least 1% trough difference to baseline
xlim([0 4]);
nhline(0, '--k');
troughTime = nan(nrAreas, length(useIdx));
for iAreas = 1 : nrAreas
    cData = postStimTime(iAreas,useIdx);
    h(iAreas) = Violin(cData, iAreas, 'ViolinColor', earlyColor, 'EdgeColor', earlyColor, 'BoxColor', [0 0 0]);
    troughTime(iAreas,useIdx) = cData;
    
    % give feedback
    disp('==================')
    fprintf('Trough time for Area %i:\n', iAreas)
    fprintf('Trough time = %.2f %c %.2f seconds \n'  , mean(cData), char(177), sem(cData));
    disp('==================')
end
ax = h(1).ViolinPlot.Parent;
axis square
ax.XTick = 1 : nrAreas;
xlabel('Area IDs')
ylabel('time from stimulus (s)');
title('Trough time');
subtitle(sprintf('Times: %.4f, %.4f, %.4f', median(postStimTime(iAreas,useIdx), 2)'));

% recovery time
nexttile
xlim([0 4]);
nhline(0, '--k');
recTime = nan(nrAreas, length(useIdx));
for iAreas = 1 : nrAreas
    cData = zeros(1, length(useIdx));
    cData(useIdx) = recoveryTime(iAreas,useIdx) - postStimTime(iAreas,useIdx);
    h(iAreas) = Violin(cData, iAreas, 'ViolinColor', earlyColor, 'EdgeColor', earlyColor, 'BoxColor', [0 0 0]);
    recTime(iAreas,:) = cData;
    
    % give feedback
    disp('==================')
    fprintf('Recovery time for Area %i:\n', iAreas)
    fprintf('Recovery time = %.2f %c %.2f seconds \n'  , mean(cData), char(177), sem(cData));
    disp('==================')
end
ax = h(1).ViolinPlot.Parent;
axis square
ax.XTick = 1 : nrAreas;
xlabel('Area IDs')
ylabel('time from stimulus (s)');
title('Recovery time');
subtitle(sprintf('Times: %.4f, %.4f, %.4f', nanmedian(recTime,2)));

%% mean movie for all recordings
cData = [];
batchSize = 500;
for iRecs = 1 : nrRecs
    cV = allV{iRecs};
    baseline = nanmean(cV(:, 1 : baseFrames),2);
    cV = insertColumns(cV, baseFrames+1, nan, gapWidth);
    cV = smoothCol(cV, 2, 15, 'box');
    cIdx = find(timeAxis >= plotRange(1) & timeAxis <= plotRange(2));
    
    if iRecs == 1
        cData = svdFrameReconstruct(allU{iRecs}, cV(:, cIdx) - baseline);
    else
        cData = cData + ((svdFrameReconstruct(allU{iRecs}, cV(:, cIdx) - baseline) - cData) / iRecs); %running average
    end
    loopReporter(iRecs, nrRecs, 10)
end

%% compute mean activity for each area
timeRange = [-300 300]; %time range for prestim and poststim period in seconds
[areaMask, areaLabels] = getAllenAreas('allenMask', 10, true);
nrAreas = length(areaLabels);
areaData = cell(2, nrAreas);
for iRecs = 1 : nrRecs
    cV = allV{iRecs};
    
    preStimStart = max([1, baseFrames + round(timeRange(1) .* opts2.sRate)]);
    postStimEnd = min([size(cV,2), baseFrames + round(timeRange(2) .* opts2.sRate)]);
    
    preStim = svdFrameReconstruct(allU{iRecs}, nanmean(cV(:, preStimStart : baseFrames),2));
    postStim = svdFrameReconstruct(allU{iRecs}, nanmean(cV(:, baseFrames + 1 : postStimEnd),2));
    
    for iAreas = 1 : nrAreas
        usePixels = areaMask{iAreas};
        usePixels(:, 1:586/2) = false;
        
        areaData{1, iAreas}(iRecs) = nanmean(preStim(usePixels));
        areaData{2, iAreas}(iRecs) = nanmean(postStim(usePixels));
    end
end

% show primary motor cortex activity
figure
clear h
cIdx = strcmpi(areaLabels, 'MOp');
preMotor = areaData{1, cIdx};
postMotor = areaData{2, cIdx};
motorData = postMotor - preMotor;
h = Violin(motorData, 0, 'ViolinColor', earlyColor, 'EdgeColor', earlyColor, 'BoxColor', [0 0 0], 'BandWidth', 0.001);
xlim([-1 1])
axis square
ax = h.ViolinPlot.Parent;
ax.XTick = 0;
ax.XTickLabel = {'PostStim'};
ylabel('dF/F (%)');
title(sprintf('Mean activity %s, %d to %d seconds\n postStim dF/F = %.4f; pVal = %.2f', ...
    areaLabels{cIdx}, timeRange(1), timeRange(2), median(motorData), signrank(motorData)));
nhline(0, 'k--')
niceFigure;
