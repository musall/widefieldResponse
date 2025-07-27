function dataOut = smoothCol(dataIn, dim, fWidth, fType, fLength)
% function to smooth columns of a data matrix 'dataIn'.
% dataOut = smoothCol(dataIn, dim, fWidth, fType, fLength)
% fType defines the filter type which is either a box, gaussian or exponential filter.
% When using a box filter 'fWidth' is the size of the moving average, when
% using a gaussian filter 'fWidth' defines its full width half maximum.
% dim is dimension that is smooth over (default is 1st dimension).
% Default filter is a 5pt box filter when only 'dataIn' is provided.
% ----------------------------------------
% To make the code run faster when using the gaussian filter, you can also
% use a shorter filter length. Default is 100*sigma.
%
% Simon Musall
% Cold Spring Harbor, 8/15/2019

if ~exist('fWidth','var') || isempty(fWidth)
    fWidth = 5; %default filter length is 5 dpoints
end

if ~exist('fType','var') || isempty(fType)
    fType = 'box'; %default is box filter
end

if ~exist('dim','var') || isempty(dim)
    dim = 1; %default is smoothing first dimension
end

if ~ismember(fType,{'box' 'gauss' 'exp' 'invExp', 'half-gauss'})
    error('Unknown filter type. Use "box", "gauss", "exp", , " or "invExp" for fType')
end

dSize = size(dataIn);
if dim > length(dSize)
    error('Dimension to be smoothed is above the dimensionality of the input data.');
end

% permute dataset to allow smoothing over other dimensions
dimOrder = 1: length(dSize);
if dim ~= 1
    dimOrder(dim) = [];
    dimOrder = [dim dimOrder];
    dataIn = permute(dataIn, dimOrder);
end

%add buffer to first dimension to avoid edge artefacts
edgeWidth = ceil(fWidth); %make sure this is an integer
cbegin = repmat(dataIn(1,:,:),[edgeWidth ones(1,length(dSize)-1)]); %start buffer
cend = repmat(dataIn(end,:,:),[edgeWidth ones(1,length(dSize)-1)]); %end buffer
dataIn = cat(1,cbegin,dataIn,cend);

dataIn = reshape(dataIn,size(dataIn,1), []); %merge other dimensions for smoothing

if sum(~isnan((abs(dataIn(:))))) > 0 %check if there is any data available
    if ~strcmpi(class(dataIn),'double') %make sure data is double
        dataIn = double(dataIn);
    end
    
    % check if data is interupted by nans and separate into smaller parts
    % for smoothing. Otherwise the will be affected by nan-edges
    if any(isnan(dataIn(:)))
        nanIdx = find(any(isnan(dataIn),2));
        nanIdx = [0; nanIdx];
    else
        nanIdx = 0; %use all data
    end
    
    dataOut = nan(size(dataIn));
    for x = 1 : length(nanIdx)
        if x ~= length(nanIdx)
            cIdx = nanIdx(x)+1:nanIdx(x+1)-1; %data for current patch
        else
            cIdx = nanIdx(x)+1:size(dataIn,1); %data for current patch
        end
        
        if ~isempty(cIdx)
            useData = dataIn(cIdx,:);
            n = size(useData,1);
            
            if n > fWidth
                if strcmpi(fType,'box') %box filter
                    fWidth = fWidth - 1 + mod(fWidth,2); %make sure filter length is odd
                    cbegin = cumsum(useData(1:fWidth-2,:),1);
                    cbegin = bsxfun(@rdivide, cbegin(1:2:end,:), (1:2:(fWidth-2))');
                    cend = cumsum(useData(n:-1:n-fWidth+3,:),1);
                    cend = bsxfun(@rdivide, cend(end:-2:1,:), (fWidth-2:-2:1)');
                    outData = conv2(useData,ones(fWidth,1)/fWidth,'full'); %smooth trace with moving average
                    outData = [cbegin;outData(fWidth:end-fWidth+1,:);cend];
                    
                elseif strcmpi(fType,'gauss') %gaussian filter
                    fSig = fWidth./(2*sqrt(2*log(2))); %in case of gaussian smooth, convert fwhm to sigma.
                    if ~exist('fLength','var') || isempty(fLength)
                        fLength = round(fSig * 100); %length of gaussian filter
                    end
                    fLength = fLength-1+mod(fLength,2); % ensure kernel length is odd
                    kernel = exp(-linspace(-fLength / 2, fLength / 2, fLength) .^ 2 / (2 * fSig ^ 2));
                    kernel = kernel / norm(kernel); %normalize kernel
                    outData = conv2(useData,kernel','same'); %smooth trace with gaussian
                    
                elseif strcmpi(fType,'half-gauss') %half-gaussian filter
                    fSig = fWidth./(2*sqrt(2*log(2))); %in case of gaussian smooth, convert fwhm to sigma.
                    if ~exist('fLength','var') || isempty(fLength)
                        fLength = round(fSig * 100); %length of gaussian filter
                    end
                    fLength = fLength-1+mod(fLength,2); % ensure kernel length is odd
                    kernel = exp(-linspace(-fLength / 2, fLength / 2, fLength) .^ 2 / (2 * fSig ^ 2));
                    kernel(1:floor(fLength/2)) = 0; %cut left side of the kernel
                    kernel = kernel / norm(kernel); %normalize kernel
                    outData = conv2(useData,kernel','same'); %smooth trace with gaussian
                    
                elseif strcmpi(fType,'exp') %exponential filter
                    if ~exist('fLength','var') || isempty(fLength)
                        fLength = ceil(fWidth * 5);
                    end
                    kernel = exponentialFilter(fWidth, fLength);
                    kernel = [zeros(1, length(kernel)), kernel]';
                    outData = conv2(useData,kernel,'same'); %convolve trace with exponential
                    
                elseif strcmpi(fType,'invExp') %inverse exponential filter for deconvolution
                    if ~exist('fLength','var') || isempty(fLength)
                        fLength = ceil(fWidth * 5);
                    end
                    kernel = exponentialFilter(fWidth, fLength);
                    kernel = [zeros(1, length(kernel)), kernel]'; % Make it a column vector for 2D conv
                    
                    % Ensure the kernel has the same size as the data
                    dataSize = size(useData);
                    kernelPadded = zeros(dataSize);
                    ks = size(kernel);
                    kernelPadded(1:ks(1), 1:ks(2)) = kernel;
                    kernelCentered = circshift(kernelPadded, -floor(size(kernel)/2)); %make sure kernel is centered
                    
                    % Compute 2D FFT of both the convolved data and the kernel
                    OUT = fft2(useData);
                    KERNEL = fft2(kernelCentered);
                    
                    % Regularization to avoid division by small numbers
                    epsilon = 1e-8;
                    outData = ifft2(OUT ./ (KERNEL + epsilon));
                    
                    % Take the real part (imaginary part should be negligible)
                    outData = real(outData);
                end
                
                % combine patches to output array
                dataOut(cIdx,:) = outData;
                
            else
                dataOut(cIdx,:) = useData;
            end
        end
    end
    
else
    dataOut = dataIn;
end

dataOut = reshape(dataOut,[size(dataOut,1) dSize(dimOrder(2:end))]); %split dimensions again
dataOut = dataOut(edgeWidth+1:end-edgeWidth, : ,:); %remove buffers
if dim ~= 1
    dataOut = ipermute(dataOut, dimOrder);
end


function kernel = exponentialFilter(decay_time, kernel_length)
% Generate a 1D exponential filter kernel with specified decay time.
% Decay time indicates after how many samples 1/e (37%) of the intial value
% is reached

% Compute the decay factor from the decay time
alpha = exp(-1 / decay_time);

% Generate the kernel values
kernel = zeros(1, kernel_length);
for t = 1:kernel_length
    kernel(t) = alpha^(t-1);
end

% Normalize the kernel to have unit energy
kernel = kernel / sum(kernel);

% Adjust the kernel length if necessary to ensure that it sums to 1
kernel_sum = sum(kernel);
if abs(kernel_sum - 1) > eps
    diff = 1 - kernel_sum;
    kernel(end) = kernel(end) + diff;
end