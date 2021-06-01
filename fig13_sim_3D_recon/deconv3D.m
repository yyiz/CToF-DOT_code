addpath("../lib");

Jdir = "3_21_21_sim_3D_Jacobian";
Jfname = "J_coloc_nSD=32x32_voxSize=32x32x6";
Jfile = "J";

gamma = 0.25;
numTicks_recon = 6;
numTicks_truth = 7;

absMua = 0.01 .* [1 2 3 4 5 6];
intTime = 1e1;
maxPhots = 5e6 * intTime;

expBinWidth = 64; % Gate width, used when summing up multiple time bins
J_layers = 1:6; % Choose which layers in Z to use for Jacobian, here using top 6 layers

useFluor = false;
tau = 1000;

step = 149.6569;
lam1 = [15e-6 30e-7 20e-7 8e-7 2.1e-7 1.2e-7];
fistaOpts.lam2 = 4e-4;
fistaOpts.maxItr = 75;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;

threshVec = [0.2, 0.2, 0.45, 0.1, 0.2, 0.2];

%% Plotting parameters

fontSizeVal = 24;

im_border = 3; % pixels
nrows = 3; % truth, recon, thresh
sliceImColor = [255, 0, 0]; % R G B values
rgb_axis_len = 100;
border_width = 1;
overall_border = 10;
dx = -0.06; dw = 0; dh = -0.3; dy = -dh/2;
cbarFmt = "%0.1e";

SLAB_MINX = 0; SLAB_MINY = 0; SLAB_MINZ = -10;
layerDz = 30;
slab_color = [0.96, 0.96, 0.96];
slab_alpha = 0.3;
slab_edge_alpha = 0.2;

shouldSave = true;
savedir = "../../figs/fig13_sim_3D_recon";
splitTruthName = split(truthImName, '_');
savename = splitTruthName(2);
savedir = sprintf("%s/deconv_%s", savedir, savename);

plotScene = true;

%% Generate and visualize scene
addpath("../lib");

rng(1);

load(sprintf("../dat/%s/%s.mat", Jdir, truthImName));
load(sprintf("../dat/%s/%s.mat", Jdir, Jfname));

if ~(useFluor) % Load background transient
    load(sprintf("../dat/%s/tpsf_compress.mat", Jdir));
end

J = J_final; clear J_final;

vars = fieldnames(Jheaders);
for i = 1:length(vars)
    assignin('base', vars{i}, Jheaders.(vars{i}));
end

mu = reshape(absMua, 1, 1, []) .* truth;
mu_orig = mu;

layerDx = VOX_W; layerDy = VOX_L; 

mu = mu(:);

%% Apply forward model

nVox = prod([VOX_L, VOX_W, VOX_H]);
nSrcDet = prod([SRC_L, SRC_W, SENS_L, SENS_W]);
if (nSrcDet ~= size(J,1)/NBINS) % Collocated measurements
    nSrcDet = SRC_L * SRC_W;
end

% Apply fluorescence lifetime
if (useFluor)
    J = reshape(J, Jheaders.NBINS, []);
    binW = (TIME_MAX - TIME_MIN) / NBINS;
    
    [Jfluor, jTimeAx] = convTpsf(J, 'exp', binW, 'lifetime', tau, 'multiCol', true);
    
    old_NBINS = NBINS;
    NBINS = length(jTimeAx);
    Jheaders.NBINS = NBINS;
    J = Jfluor;
end

m = J*mu; % Generate ground truth measurements

m_cube = permute(reshape(m, NBINS, VOX_L, VOX_W), [2 3 1]);

f_intM_Jacobian = figure('Position', [115.5 351.5 604.5 518]); 
imagesc(sum(m_cube,3));  cb_Jmu = colorbar;
% title("$m=J\mu$", 'interpreter', 'latex');
set(gca, 'FontSize', fontSizeVal);
set(gca, 'FontName', 'Times New Roman');
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); axis image;

%%

Jfile_struct = load(sprintf("../dat/%s/%s.mat", Jdir, Jfile));

J_conv = Jfile_struct.J;
J_conv = accumTimeBin(J_conv, expBinWidth);
J_conv = reshape(J_conv, NBINS, 64, 64, 10);
J_conv = J_conv(:,:,:,J_layers);
J_psf = permute(J_conv, [2 3 1 4]);

f_fista = @(X)fwd(X, J_psf);
ft_fista = @(Y)adj(Y, J_psf);

m_conv = f_fista(mu_orig); % Apply forward model using conv approx

f_intM_conv = figure('Position', [115.5 351.5 604.5 518]); 
imagesc(sum(m_conv,3));
cb_conv = colorbar;
set(gca, 'FontSize', fontSizeVal);
set(gca, 'FontName', 'Times New Roman');
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); axis image;

%% Apply noise to measurements

m_conv_reshape = transpose(reshape(m_conv, [], NBINS));
intBkgSig = sum(bkgTpsf(:));
scaleFac = maxPhots ./ intBkgSig;

bkgMeas = repmat(bkgTpsf, [1, nSrcDet]);
absMeas = abs(bkgMeas - m_conv_reshape); % Handle numerical imprecision with abs

bkgMeas = bkgMeas .* scaleFac;
absMeas = absMeas .* scaleFac;
m_noisy = poissrnd(bkgMeas) - poissrnd(absMeas);

m_noisy = m_noisy ./ scaleFac;

m_noisy = permute(reshape(m_noisy, NBINS, SRC_L, SRC_W), [2 3 1]);

%% Test image reconstruction

sizeX = [VOX_L VOX_W VOX_H];

lam1Vec = ones(VOX_L, VOX_W, VOX_H);
lam1Vec = reshape(lam1, 1, 1, []) .* lam1Vec;

fprintf("Running FISTA...");
tic; [fistaRecon, ~] = fista_layerReg(m_noisy,f_fista,ft_fista,step,sizeX,fistaOpts,lam1Vec); fistaRuntime = toc;
fprintf("done! Finished in: %.3f seconds\n", fistaRuntime);

% Apply threshold to reconstructed images, using threshVal
thresh = fistaRecon;
for i = 1:size(thresh,3)
    layer_i = thresh(:,:,i);
    layer_i = layer_i ./ max(layer_i(:));
    threshInds = layer_i > threshVec(i);
    layer_i(threshInds) = 1;
    layer_i(~threshInds) = 0;
    thresh(:,:,i) = layer_i;
end

%% Plot results

% Plot scene
ftruth = plotRecon3D(mu_orig, 'rgb', sliceImColor, 'cbar', false, 'normalizeIm', false);
view([227 35]);
cb = colorbar;
maxPix = max(mu_orig(:));
cb.Ticks = linspace(0,maxPix,numTicks_truth);
tlabels = [];
for i = 1:numTicks_truth
    tlabels = [tlabels sprintf(cbarFmt, cb.Ticks(i))];
end
cb.TickLabels = tlabels;
pos =  cb.Position; %gets the positon and size of the color bar
set(cb,'Position',[pos(1)+dx pos(2)+dy pos(3)+dw pos(4)+dh])% To change size
set(gca, 'FontSize', fontSizeVal);
set(gca, 'FontName', 'Times New Roman');
set(gca,'visible', 'off');
voxel([SLAB_MINX, SLAB_MINY, SLAB_MINZ],[layerDx, layerDy, layerDz],slab_color,slab_alpha, slab_edge_alpha);


% Plot reconstruction
fistaRecon_hdr = fistaRecon .^ gamma;
frecon = plotRecon3D(fistaRecon_hdr, 'rgb', sliceImColor, 'cbar', false, 'normalizeIm', false);
view([227 35]);
cb = colorbar;
maxPix = max(fistaRecon_hdr(:));
cb.Ticks = linspace(0,maxPix,numTicks_recon);
tlabels = [];
for i = 1:numTicks_recon
    tlabels = [tlabels sprintf(cbarFmt, cb.Ticks(i) .^ (1/gamma))];
end
cb.TickLabels = tlabels;
pos =  cb.Position; %gets the positon and size of the color bar
set(cb,'Position',[pos(1)+dx pos(2)+dy pos(3)+dw pos(4)+dh])% To change size
set(gca, 'FontSize', fontSizeVal);
set(gca, 'FontName', 'Times New Roman');
set(gca,'visible', 'off');
voxel([SLAB_MINX, SLAB_MINY, SLAB_MINZ],[layerDx, layerDy, layerDz],slab_color,slab_alpha, slab_edge_alpha);


% Plot threshold image
fthresh = plotRecon3D(thresh, 'rgb', sliceImColor);
view([227 35]);
set(gca, 'FontSize', fontSizeVal);
set(gca, 'FontName', 'Times New Roman');
xlabel(''); ylabel(''); zlabel('');
colorbar('off');
set(gca,'visible', 'off');
voxel([SLAB_MINX, SLAB_MINY, SLAB_MINZ],[layerDx, layerDy, layerDz],slab_color,slab_alpha, slab_edge_alpha);


% Plot Z-slices
plotLayers = 1:6;

imH = 2*overall_border + im_border*(length(plotLayers)+1) + VOX_W*length(plotLayers);
imW = 2*overall_border + im_border*(nrows+1) + nrows*VOX_L;
slicesIm = ones(imH, imW, 3); % 3rd dimension for RGB

rax = linspace(1, sliceImColor(1)./255, rgb_axis_len)';
gax = linspace(1, sliceImColor(2)./255, rgb_axis_len)';
bax = linspace(1, sliceImColor(3)./255, rgb_axis_len)';

sliceImColor = reshape(sliceImColor, 1, 1, []);
for i = 1:length(plotLayers)
    truthCol = (0)*(im_border+VOX_W) + im_border + overall_border;
    truthRow = (i-1)*(im_border+VOX_W) + im_border + overall_border;
    truthRowEnd = truthRow+VOX_L-1;
    truthColEnd = truthCol+VOX_W-1;
    % Add black border
    slicesIm((truthRow-border_width):(truthRowEnd+border_width),...
        (truthCol-border_width):(truthColEnd+border_width),:) = 0;
    
    truthImTemp = mu_orig(:,:,plotLayers(i)); 
    truthImTemp = truthImTemp ./ max(truthImTemp(:));
    rslice = interp1(linspace(0, 1, rgb_axis_len), rax, truthImTemp);
    gslice = interp1(linspace(0, 1, rgb_axis_len), gax, truthImTemp);
    bslice = interp1(linspace(0, 1, rgb_axis_len), bax, truthImTemp);
    truthIm = cat(3, rslice, gslice, bslice);
    slicesIm(truthRow:truthRowEnd, truthCol:truthColEnd, 1:3) = truthIm;


    reconCol = (1)*(im_border+VOX_W) + im_border + overall_border;
    reconRow = (i-1)*(im_border+VOX_W) + im_border + overall_border;
    reconRowEnd = reconRow+VOX_L-1;
    reconColEnd = reconCol+VOX_W-1;
    % Add black border
    slicesIm((reconRow-border_width):(reconRowEnd+border_width),...
        (reconCol-border_width):(reconColEnd+border_width),:) = 0;
    
    reconImTemp = fistaRecon(:,:,plotLayers(i));  
    reconImTemp = reconImTemp ./ max(reconImTemp(:));
    rslice = interp1(linspace(0, 1, rgb_axis_len), rax, reconImTemp);
    gslice = interp1(linspace(0, 1, rgb_axis_len), gax, reconImTemp);
    bslice = interp1(linspace(0, 1, rgb_axis_len), bax, reconImTemp);
    reconIm = cat(3, rslice, gslice, bslice);
    slicesIm(reconRow:reconRowEnd, reconCol:reconColEnd, 1:3) = reconIm;

    
    threshCol = (2)*(im_border+VOX_W) + im_border + overall_border;
    threshRow = (i-1)*(im_border+VOX_W) + im_border + overall_border;
    threshRowEnd = threshRow+VOX_L-1;
    threshColEnd = threshCol+VOX_W-1;
    % Add black border
    slicesIm((threshRow-border_width):(threshRowEnd+border_width),...
        (threshCol-border_width):(threshColEnd+border_width),:) = 0;

    threshImTemp = thresh(:,:,plotLayers(i));
    threshImTemp = threshImTemp ./ max(threshImTemp(:));
    rslice = interp1(linspace(0, 1, rgb_axis_len), rax, threshImTemp);
    gslice = interp1(linspace(0, 1, rgb_axis_len), gax, threshImTemp);
    bslice = interp1(linspace(0, 1, rgb_axis_len), bax, threshImTemp);
    threshIm = cat(3, rslice, gslice, bslice);
    slicesIm(threshRow:threshRowEnd, threshCol:threshColEnd, 1:3) = threshIm;

    fprintf("Depth: %f\n", VOX_ORIG(3) + (plotLayers(i)-1)*Jheaders.VOX_ZLEN);
end
fprintf("Full range: %f-%f\n", Jheaders.VOX_ORIGZ, Jheaders.VOX_ORIGZ + (Jheaders.VOX_H * Jheaders.VOX_ZLEN));

fslices = figure('Position', [155.5 107.5 1184 815.5]); imagesc(slicesIm);
axis image
set(gca, 'visible', 'off');

%% Save results

if (shouldSave)
    
if ~(exist(savedir, 'dir'))
    mkdir(savedir);
end
    
truthFigPath = sprintf("%s/truth_%s.fig", savedir, savename);
truthPngPath = sprintf("%s/truth_%s.png", savedir, savename);

reconFigPath = sprintf("%s/recon_%s.fig", savedir, savename);
reconPngPath = sprintf("%s/recon_%s.png", savedir, savename);

threshFigPath = sprintf("%s/thresh%s.fig", savedir, savename);
threshPngPath = sprintf("%s/thresh_%s.png", savedir, savename);

slicesFigPath = sprintf("%s/slices_%s.fig", savedir, savename);
slicesPngPath = sprintf("%s/slices_%s.png", savedir, savename);

intM_J_figPath = sprintf("%s/intM_J_%s.fig", savedir, savename);
intM_J_pngPath = sprintf("%s/intM_J_%s.png", savedir, savename);

intM_conv_figPath = sprintf("%s/intM_conv_%s.fig", savedir, savename);
intM_conv_pngPath = sprintf("%s/intM_conv_%s.png", savedir, savename);


savefig(ftruth, truthFigPath);
export_fig(ftruth, truthPngPath, '-m3', '-transparent', '-png');

savefig(frecon, reconFigPath);
export_fig(frecon, reconPngPath, '-m3', '-transparent', '-png');

savefig(fthresh, threshFigPath);
export_fig(fthresh, threshPngPath, '-m3', '-transparent', '-png');

savefig(fslices, slicesFigPath);
export_fig(fslices, slicesPngPath, '-m3', '-transparent', '-png');

savefig(f_intM_Jacobian, intM_J_figPath);
export_fig(f_intM_Jacobian, intM_J_pngPath, '-m3', '-transparent', '-png');

savefig(f_intM_conv, intM_conv_figPath);
export_fig(f_intM_conv, intM_conv_pngPath, '-m3', '-transparent', '-png');

end

%% Forward and adjoint functions

function fwdout = fwd(fwdin, psf)
    NBINS = size(psf, 3);
    VOX_L = size(fwdin, 1);
    VOX_W = size(fwdin, 2);
    VOX_H = size(fwdin, 3);
    
    fwdin_pad = padarray(fwdin, [VOX_L/2, VOX_W/2]);
    
    fwdout = zeros(VOX_L, VOX_W, NBINS);
    
    rowInd_start = VOX_L/2 + 1;
    rowInd_end = rowInd_start + VOX_L - 1;
    colInd_start = VOX_W/2 + 1;
    colInd_end = colInd_start + VOX_W - 1;
    
    for i = 1:VOX_H
        psf_layeri = psf(:,:,:,i);
        fwd_layeri = fwdin_pad(:,:,i);
        
        filt_layeri = fftFilt2(fwd_layeri, psf_layeri);
        
        fwdout = fwdout + filt_layeri(rowInd_start:rowInd_end,colInd_start:colInd_end,:);
    end
end

function adjout = adj(adjin, psf)
    VOX_H = size(psf, 4);
    VOX_L = size(adjin, 1);
    VOX_W = size(adjin, 2);
    
    adjin_pad = padarray(adjin, [VOX_L/2, VOX_W/2]);
    
    adjout = zeros(VOX_L, VOX_W, VOX_H);
    
    rowInd_start = VOX_L/2 + 1;
    rowInd_end = rowInd_start + VOX_L - 1;
    colInd_start = VOX_W/2 + 1;
    colInd_end = colInd_start + VOX_W - 1;
    
    for i = 1:VOX_H
        psf_layeri = psf(:,:,:,i);
        filt_layeri = sum(fftFilt2(adjin_pad, rot90(psf_layeri, 2)), 3);        
        adjout(:,:,i) = filt_layeri(rowInd_start:rowInd_end,colInd_start:colInd_end);
    end    
end
