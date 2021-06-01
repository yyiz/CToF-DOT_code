% MULTIPLEXING EXPERIMENTAL DATA
sizeDat = size(bkg_data_temp); sizeDat = sizeDat(1:5);
n = size(bkg_data_temp,1);
N = n*n;
H = hadamard(N);
HT = H';
bkg_multiplex = zeros(sizeDat);
meas_multiplex = zeros(sizeDat);
for i = 1:N

    Hi = reshape(H(i,:), [n, n]);
    [srcRow, srcCol] = ind2sub([n, n], i);

    bkg_i = sum(sum(Hi .* bkg_data_temp(:,:,:,:,:,i), 1), 2);
    bkg_multiplex(srcCol,srcRow,:,:,:) = bkg_i;

    meas_i = sum(sum(Hi .* measure_data_temp(:,:,:,:,:,i), 1), 2);
    meas_multiplex(srcCol,srcRow,:,:,:) = meas_i;
end
% DEMULTIPLEX
bkg_data_temp = zeros(sizeDat(1:2));
measure_data_temp = zeros(sizeDat(1:2));
for i = 1:N
    HTi = reshape(HT(i,:), [n, n]);
    [srcRow, srcCol] = ind2sub([n, n], i);

    bkg_data_temp(srcCol, srcRow,:,:,:) = sum(HTi .* squeeze(bkg_multiplex(:,:,srcCol, srcRow)), [1, 2]);
    measure_data_temp(srcCol, srcRow,:,:,:) = sum(HTi .* squeeze(meas_multiplex(:,:,srcCol, srcRow)), [1, 2]);
end