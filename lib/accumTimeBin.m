% Assume time along dimension 1

function matOut = accumTimeBin(matIn, bin_fac)

nrows = size(matIn, 1);
nrowsOut = ceil(nrows/bin_fac);

numRowsAppend = bin_fac - mod(nrows, bin_fac);
if (numRowsAppend == bin_fac)
    numRowsAppend = 0;
end
matOut = [matIn; zeros(numRowsAppend, size(matIn, 2))];

matOut = sum(reshape(matOut, bin_fac, []), 1);
matOut = reshape(matOut, nrowsOut, []);

end