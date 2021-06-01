function [f, ft, step] = mat2Handle(A, type, varargin)

    name_list = {'VOX_L', 'VOX_W', 'lvl', 'nFilts'};
    optargs = {0, 0, 0, 1};
    optargs = parseVarargin(varargin, name_list, optargs);
    [VOX_L, VOX_W, lvl, nFilts] = optargs{:};

    if (startsWith(type, 'mvp')) % matrix vector product
        if (strcmp(type, "mvp")) % 'mvp'
            AT = transpose(A);
            f = @(x) A*x;
            ft = @(y) AT*y;
        elseif (endsWith(type, "dct")) % 'mvp-dct'
            AT = transpose(A);
            f = @(x) dctBlur(A, x);
            ft = @(y) dctBlurAdj(y, AT, VOX_L, VOX_W);
        elseif (endsWith(type, "dwt")) % 'mvp-dwt'
            AT = transpose(A);
            f = @(x) waveletBlur(A, x, lvl);
            ft = @(y) waveletBlurAdj(y, AT, VOX_L, VOX_W, lvl);
        else
            error("Invalid function handle type");
        end
        
        step = 1/(2*eigs(A'*A, 1));
        
    elseif (startsWith(type, "conv")) % A is the psf in this case
        
        f = @(I) fftFilt2(repmat(I, [1, 1, nFilts]), A);
        ft = @(I) sum(fftFilt2(I, rot90(A,2)), 3);

        steps = zeros(size(A,3), 1);
        for i = 1:size(A, 3)
            cRow = ceil(size(A(:,:,i), 1)/2); cCol = ceil(size(A(:,:,i), 2)/2);
            Ieigs = fft2(circshift(A(:,:,i), 1 - [cRow, cCol]));
            Ieigs2 = Ieigs.^2;
            steps(i) = 1./(2*max(Ieigs2(:))); % step = 1/(2*lam_max(A'*A))
        end
        step = min(steps);
        
    else
        error("Invalid function handle type");
    end
end

function y = dctBlur(A, x)
    xVec = idct2(x);
    xVec = reshape(xVec, numel(xVec), 1);
    y = A*xVec;
end

function x = dctBlurAdj(y, AT, VOX_L, VOX_W)
    xVec = AT*y;
    x = dct2(reshape(xVec, VOX_L, VOX_W));
end

function y = waveletBlur(A, x, lvl)
    xIm = wavelet('2D Haar', -lvl, x);
    y = A*xIm(:);
end

function x = waveletBlurAdj(y, AT, VOX_L, VOX_W, lvl)
    xVec = AT*y;
    xIm = reshape(xVec, [VOX_L VOX_W]);
    x = wavelet('2D Haar', lvl, xIm);
end