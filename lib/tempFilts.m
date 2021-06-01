function [J_mc, m] = tempFilts(J, diffTpsf, Jheaders, filts, varargin)
    
    name_list = {'s', 'filtWeights', 'fwdModelType'};
    optargs = {0.001, [1 1 1], 'mvp'};
    optargs = parseVarargin(varargin, name_list, optargs);
    [s, filtWeights, fwdModelType] = optargs{:};

    % Temporal filters for Jacobians
    nSrcDet = prod([Jheaders.SRC_L, Jheaders.SRC_W, Jheaders.SENS_L, Jheaders.SENS_W]);
    nVox = prod([Jheaders.VOX_L, Jheaders.VOX_W, Jheaders.VOX_H]);
    
    if (isempty(J))
        J_mc = [];
    else
        jMat = reshape(J, [Jheaders.NBINS, nSrcDet, nVox]);

        timeAx = linspace(Jheaders.TIME_MIN, Jheaders.TIME_MAX, Jheaders.NBINS+1);
        timeAx = transpose(timeAx(1:end-1));

        J_mc = [];
        if contains(filts, "E")
            Eind = strfind(filts, 'E');
            Eweight = filtWeights(Eind);
            Efilt_J = Eweight*squeeze(trapz(timeAx, jMat, 1));
            J_mc = [J_mc; Efilt_J];
        end
        if contains(filts, "M")
            Mind = strfind(filts, 'M');
            Mweight = filtWeights(Mind);
            Mfilt_J = Mweight*squeeze(trapz(timeAx, timeAx .* jMat, 1));
            J_mc = [J_mc; Mfilt_J];
        end
        if contains(filts, "L")
            Lind = strfind(filts, 'L');
            Lweight = filtWeights(Lind);
            Lfilt_J = Lweight*squeeze(trapz(timeAx, exp(-s.*timeAx).*jMat, 1));
            J_mc = [J_mc; Lfilt_J];
        end
    end
        
    m = [];
    if (isempty(diffTpsf))
        
    else
        
        if (strcmp(fwdModelType, 'conv')) % if using a convolutional model
            nSrcDet = nVox;
        end
        
        if (isfield(Jheaders, 'nMeasBins')) % Temporal filters for experimental measurements
            diffTpsf = reshape(diffTpsf, [Jheaders.nMeasBins, nSrcDet]);
            measTimeAx = Jheaders.measTimeAx;
        else
            diffTpsf = reshape(diffTpsf, [Jheaders.NBINS, nSrcDet]);
            measTimeAx = timeAx;
        end
        
        if contains(filts, "E")
            Efilt_diff = transpose(trapz(measTimeAx, diffTpsf, 1));
            m = [m; Efilt_diff];
        end
        if contains(filts, "M")
            Mfilt_diff = transpose(trapz(measTimeAx, measTimeAx .* diffTpsf, 1));
            m = [m; Mfilt_diff];
        end
        if contains(filts, "L")
            Lfilt_diff = transpose(trapz(measTimeAx, exp(-s.*measTimeAx).*diffTpsf, 1));
            m = [m; Lfilt_diff];
        end
    end
    
end