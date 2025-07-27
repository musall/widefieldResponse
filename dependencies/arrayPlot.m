function cLine = arrayPlot(amatrix,varargin)
% plot all columns of a data matrix as individual lines and add the mean
% over columns as a single thick line. aside of the data matrix, additional
% inputs that are compatible with the 'plot' command can be provided.
% usage: cLine = arrayPlot(amatrix,varargin)

%check if first input is for x-axis
if ~isempty(varargin) && numel(amatrix) == size(varargin{1},1)
    xAxis = amatrix;
    amatrix = varargin{1};
    varargin = varargin(2:end);
else
    xAxis = 1 : size(amatrix,1);
end

% plot individual columns
cLine = plot(xAxis, amatrix, varargin{:});

% check if hold is already on
a = gca;
if ishold(a)
    checker = true;
else
    hold(a,'on'); checker = false;
end

% get color and linewidth
cLineWidth = cLine(1).LineWidth;
cColor = cLine(1).Color + 0.75;
cColor(cColor >  1) = 1;
for x = 1 : length(cLine)
    cLine(x).Color = cColor;
end

% show mean over all columns
cLine = plot(xAxis, nanmean(amatrix,2), varargin{:}, 'linewidth', cLineWidth * 4);
  
if ~checker
    hold(a,'off');
end


