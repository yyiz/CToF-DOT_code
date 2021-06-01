rng(1);

plotTacq = [1, 5, 15, 20];

load("../dat/11_11_20_multiplex/J5.mat");
truthFname = "R2";
bw_thresh = 60; % When converting png to binary image, this is pixel thresholding (0-255)
dcr = 200;

savefname = "bkgnoise_psnr_vs_compression_multiplexing";

N = Jheaders.SRC_L;
nbins = Jheaders.NBINS;
nvox = prod([Jheaders.VOX_L, Jheaders.VOX_W, Jheaders.VOX_H]);
arrDim = N*N;
H = hadamard(arrDim);

fwdType = 'mvp';
fistaOpts.lam1 = 0;
fistaOpts.lam2 = 0;
fistaOpts.maxItr = 100;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;

intTimeVec = linspace(1e1,1e4,30);

%% Generate scene and simulated measurements

addpath("../lib");
addpath("../export_fig");


defaultPath = "../dat";
savedir = "../../figs/fig10_multiplexing_results";

if ~exist(savedir, 'dir')
    mkdir(savedir);
end

delUa = 0.01;
truthIm = sprintf("truth_ims/%s.png", truthFname);
mu = imread(sprintf("%s/%s", defaultPath, truthIm));
mu = imresize(mu, [Jheaders.VOX_L, Jheaders.VOX_W]);
if (size(mu,3) > 1)
    mu = rgb2gray(mu);
end
mu = double(mu);
absInds = mu < bw_thresh;
mu(absInds) = delUa;
mu(~absInds) = 0;
truthIm = double(mu);
truthIm = 255 * (truthIm ./ max(truthIm(:)));
mu = double(mu(:));
m = J*mu;

% Restructure data to (sc, sr, dc, dr,nbins)
mfull = mat2dat(m, Jheaders, false);
Jfull = mat2dat(J, Jheaders, true);
bkgfull = mat2dat(bkgTpsf, Jheaders, false);

% Handle numerical imprecision
mfull(abs(mfull) < 1e-18) = 0;

Jfull_coloc = zeros(N, N, nbins, nvox);
mfull_coloc = zeros(N, N, nbins);
bkg_coloc = zeros(N, N, nbins);
for nr = 1:N
    for nc = 1:N
        Jfull_coloc(nc, nr, :, :) = Jfull(nc, nr, nc, nr, :, :);
        mfull_coloc(nc, nr, :) = mfull(nc, nr, nc, nr, :);
        bkg_coloc(nc, nr, :) = bkgfull(nc, nr, nc, nr, :);
    end
end

no_multiplex_psnr_vec = zeros(1,length(intTimeVec));

for i = 1:length(intTimeVec)
    
intTime = intTimeVec(i);
maxPhots = 5e3 * intTime; % Assume pile-up limit is at 5000 counts per millisecond
    
% Apply poisson noise
int_bkg_coloc = sum(bkg_coloc, 3);
max_sig_coloc = max(int_bkg_coloc(:));
scaleFac = maxPhots ./ max_sig_coloc;

absMeas = bkg_coloc - mfull_coloc;
bkgMeas = bkg_coloc .* scaleFac;
absMeas = absMeas .* scaleFac;
mfull_noisy = (poissrnd(bkgMeas) - poissrnd(absMeas));

mfull_noisy = mfull_noisy ./ scaleFac;

% Reshape J and m to p x q and p x 1, respectively
J_coloc = reshape(Jfull_coloc, N*N*nbins, nvox);
m_coloc = mfull_noisy(:);

% Perform image reconstruction
[f_fista, fT_fista, fistaStep] = mat2Handle(J_coloc, fwdType,...
    'VOX_L', Jheaders.VOX_L, 'VOX_W', Jheaders.VOX_W);

fprintf("Running FISTA...");
tic; [fistaRecon, fistaErr] = fista(m_coloc,f_fista,fT_fista,fistaStep,[nvox,1],fistaOpts); fistaRuntime = toc;
fprintf("done! Finished in: %.3f seconds\n", fistaRuntime);

% Don't add to PSNR plot if resulting image is all zeros
fistaRecon = reshape(fistaRecon, Jheaders.VOX_L, Jheaders.VOX_W);
if (sum(fistaRecon, 'all') > 0)
    fistaRecon = 255*(fistaRecon ./ max(fistaRecon(:)));
else
    
end

psnrVal = psnr(uint8(fistaRecon), uint8(truthIm));
no_multiplex_psnr_vec(i) = psnrVal;

if ismember(i, plotTacq)
    reconF = figure();
    imagesc(fistaRecon);
    pbaspect([1 1 1]);
    set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
    colormap(flipud(gray));
    
    fname = strrep(sprintf("%dNO_multiplex_recon_Tacq=%0.2fms", find(plotTacq==i), intTimeVec(i)), '.', 'p');
    reconName = sprintf("%s/%s.png", savedir, fname);
end

end
%% Test multiplexing

multiplex_psnr_vec = zeros(1,length(intTimeVec));

for j = 1:length(intTimeVec)

intTime = intTimeVec(j);
maxPhots = 5e3 * intTime; % Assume pile-up limit is at 5000 counts per millisecond
    
% Perform multiplexing
Hm = zeros(size(mfull));
for i = 1:arrDim
    Hi = reshape(H(i,:), N, N);
    Hi_pos = Hi > 0;
    Hi_neg = Hi < 0;

    [sc, sr] = ind2sub([N, N], i);

    bkg_pos = squeeze(sum(Hi_pos .* sum(bkgfull, 5), [1, 2]));
    bkg_neg = squeeze(sum(Hi_neg .* sum(bkgfull, 5), [1, 2]));
    max_sig_full = max([bkg_pos(:); bkg_neg(:)]);
    
    scaleFac = maxPhots ./ max_sig_full;
    
    absMeas = abs(bkgfull - mfull);
    bkgscale = bkgfull .* scaleFac;
    absscale = absMeas .* scaleFac;
    Hm_pos = poissrnd(Hi_pos .* bkgscale) - poissrnd(Hi_pos .* absscale);
    Hm_neg = poissrnd(Hi_neg .* bkgscale) - poissrnd(Hi_neg .* absscale);
    
    noisy_m = (Hm_pos - Hm_neg) ./ scaleFac;
   
    Hm(sc, sr, :, :, :) = sum(noisy_m, [1, 2]);
end

% Demultiplex (retrieve only confocal data)
HTHm = zeros(N, N, nbins);
for i = 1:arrDim
    HTi = reshape(H(:,i), N, N);
    [sc, sr] = ind2sub([N, N], i);
    % Note: by indexing a specific detector and all detectors, this is
    % de-multiplexing the sources, i.e. the HTi is being multiplied along
    % the SOURCE dimensions
    HTHm(sc, sr, :) = sum(sum(HTi .* Hm(:,:, sc, sr, :), 1), 2);
end
HTHm = HTHm ./ arrDim;
HTHm = HTHm(:);

fprintf("Running FISTA...");
tic; [fistaReconMulti, fistaErr] = fista(HTHm,f_fista,fT_fista,fistaStep,[nvox,1],fistaOpts); fistaRuntime = toc;
fprintf("done! Finished in: %.3f seconds\n", fistaRuntime);

fistaReconMulti = reshape(fistaReconMulti, Jheaders.VOX_L, Jheaders.VOX_W);

if (sum(fistaReconMulti, 'all') > 0)
    fistaReconMulti = 255*(fistaReconMulti ./ max(fistaReconMulti(:)));
else
    
end

psnrVal = psnr(uint8(fistaReconMulti), uint8(truthIm));
multiplex_psnr_vec(j) = psnrVal;

if ismember(j, plotTacq)
    reconF = figure();
    imagesc(fistaReconMulti);
    pbaspect([1 1 1]);
    set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
    colormap(flipud(gray));

    fname = strrep(sprintf("%dmultiplex_recon_Tacq=%0.2fms", find(plotTacq==j), intTimeVec(j)), '.', 'p');
    reconName = sprintf("%s/%s.png", savedir, fname);
    
end

end

%% Plot result

f1 = figure('Position', [527 499.5 625.5 210.5]);

scatter(intTimeVec, no_multiplex_psnr_vec, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(intTimeVec, multiplex_psnr_vec, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.5);

fit_no_multiplex = fit(intTimeVec', no_multiplex_psnr_vec', 'exp2');
h = plot(fit_no_multiplex, 'r');
set(h, 'LineWidth', 2);

fit_multiplex = fit(intTimeVec', multiplex_psnr_vec', 'exp2');
h2 = plot(fit_multiplex, 'b');
set(h2, 'LineWidth', 2);

L = legend();
set(L,'visible','off')

xlabel("Integration time (ms)");
ylabel("PSNR");
set(gca, 'XDir','reverse')
set(gca, 'FontSize', 14);
set(gca, 'FontName', 'Times New Roman');

%% Save result

saveName = sprintf("%s/%s.png", savedir, savefname);

if exist(saveName, 'file')
    fprintf("file already exists. Overwrite?\n");
    keyboard;
end
export_fig(f1, saveName, '-transparent', '-m3', '-png');

f2 = figure();
imagesc(truthIm);
pbaspect([1 1 1]);
set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
colormap(flipud(gray));

