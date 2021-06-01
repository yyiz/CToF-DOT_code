function [mu, crow, ccol] = genLinesScene(sep, VOX_L, VOX_W, absCoeff, varargin)

    name_list = {'lineLen', 'nreps'};
    optargs = {15, 3};
    optargs = parseVarargin(varargin, name_list, optargs);
    [lineLen, nreps] = optargs{:};
    mu = zeros(VOX_L, VOX_W);
    
    crow = round((VOX_L-lineLen)/2);
    ccol = round(VOX_W/2)-nreps*sep;
    for i = 1:nreps
        colStart = ccol + (i-1)*2*sep;
        colEnd = colStart + sep - 1;
        mu(crow:(crow+lineLen-1),colStart:colEnd)=1;
        mu(crow:(crow+lineLen-1),(colStart:colEnd)+sep)=0;
    end
    mu = mu*absCoeff;
end