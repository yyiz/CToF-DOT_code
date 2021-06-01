function [tpsfOut, tpsfAx] = convTpsf(tpsfIn, psfType, binW, varargin)

name_list = {'lifetime', 'duration', 'stdDev', 'multiCol'};
optargs = {1, 1, 1, false};
optargs = parseVarargin(varargin, name_list, optargs);
[lifetime, duration, stdDev, multiCol] = optargs{:};

nBins = size(tpsfIn, 1);

if (strcmp(psfType,'exp'))
    
    nPtsPsf = nBins * 15 * (lifetime/1000);
    maxTimeExp = binW * (nPtsPsf-1);
    timeAxExp = transpose(linspace(0, maxTimeExp, nPtsPsf));
    psfDistrib = (1/lifetime)*exp(-timeAxExp/lifetime);
    
elseif (strcmp(psfType,'gauss'))
    
    nPtsPsf = ceil((6*stdDev) ./ binW);
    psfDistrib = normpdf(transpose(1:nPtsPsf), nPtsPsf/2, nPtsPsf/6);
    
elseif (strcmp(psfType,'box'))
    
    nPtsPsf = ceil(duration / binW);
    psfDistrib = ones(nPtsPsf, 1) ./ nPtsPsf;
    
else
    error("Incorrect PSF specification");
end

psfDistrib = psfDistrib / sum(psfDistrib(:)); %Normalize expDist to 1 for convolution

if (multiCol)
    tpsfOut = conv2(psfDistrib, 1, tpsfIn);
else
    tpsfOut = conv(tpsfIn, psfDistrib);
end
tpsfAx = transpose(linspace(0, binW * (nPtsPsf+nBins-2), nPtsPsf+nBins-1));

end