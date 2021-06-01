clear; close all; clc;

%% Set variables

addpath('../lib');
addpath('../export_fig');

defaultPath = '../dat';

truthPath = "../dat/3_9_21_exp_3D_recon_data_absorption/truth/truth.mat";

loadJdir = '3_12_21_recon3D-64x64x10vox-1srcdet';
loadJname = 'J';
loadBkgName = 'tpsf';
loadJpath = sprintf("%s/%s/%s.mat", defaultPath, loadJdir, loadJname);
loadBkgPath = sprintf("%s/%s/%s.mat", defaultPath, loadJdir, loadBkgName);

expFname = "recon3D_3pillars_coloc_09-Mar-2021_17-17-32";
expdir = "3_9_21_exp_3D_recon_data_absorption";

shouldSave = true;
savedir = "../../figs/fig14_exp_3D_recon";

% Experimental data processing parameters
normalizeTransients = false; % Toggle whether normalize bkg/abs TPSF to have same max
maxk_norm = 20; % When normalizing noisy data, better to consider median of max-k points as maximum
timeWindow = 2001:6000; % Crop window for experimental TPSF, for removing zeros, not gating transients
alignInds = -20:0.1:20; % Shift points when aligning bkg/meas TPSF
expBinWidth = 200; % Number of time bins that are integrated into 1 bin (for exp data)
alignBkgInds = [16, 16];

% IRF calculation parameters
imgaussSig = 10; % Standard deviation for smoothing transient when determining IRF
bkgDeconvInds = [16,16]; % Index of experimental TPSF for deconvolution
irf_retrieval_fistaOpts.lam1 = 0; % Fista parameters for deconvolution
irf_retrieval_fistaOpts.lam2 = 0;
irf_retrieval_fistaOpts.maxItr = 10;
irf_retrieval_fistaOpts.tol = 0;
irf_retrieval_fistaOpts.nonneg = true; 
irf_retrieval_fistaOpts.showFigs = false;


% Parameters for generating collocated Jacobian from single source-detector pair:
gamma = 1.0; % Apply exponential biasing
nSD_L = 32; % Number of rows of in sd array
nSD_W = 32; % Number of columns in sd array
sdArr_Len = 32; % Physical length covered by sd array
sdArr_Width = 32; % Physical width covered by sd array
voxSize_L = 32; % Length of the voxel array (cropped from Jacobian)
voxSize_W = 32; % Width of the voxel array (cropped from Jacobian)

% Experimental image reconstruction
fwd_model = 'mvp';
fistaOpts.lam2 = 0;
fistaOpts.maxItr = 10;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;
lam1 = [1e10, 1e10, 2.5e8, 1.3e8, 4e7, 3e7, 9.5e6, 6e6, 3e6, 1.5e6];
threshVec = [0.5, 0.7, 0.5, 0.5, 0.3, 0.3, 0.16, 0.16, 0.2, 0.2];

%% Plotting parameters

% Toggle whether to plot
plotTestIrf = false;
visProcessedJ = false;
plotFistaErr = false;

recon_gamma = 0.3; % Gamma correction on reconstructed image

fontSizeVal = 19;
plotLayers = [3, 4, 5, 6, 7, 8];

% Slice plotting parameters
im_border = 3; % pixels
nPlotRows = 3; % truth, recon, thresh
sliceImColor = [255, 0, 0]; % R G B values
rgb_axis_len = 100;
border_width = 1;
overall_border = 10;

SLAB_MINX = 0; SLAB_MINY = 0; SLAB_MINZ = -10;
layerDz = 30;
slab_color = [0.96, 0.96, 0.96];
slab_alpha = 0.3;
slab_edge_alpha = 0.2;

%% Clean up experimental data: noramlize, align, and compress transients

if ~exist(savedir, 'dir')
    mkdir(savedir);
end

savenameList = split(expFname, '_'); savenameBase = savenameList(2);

loadExpPath = sprintf("%s/%s/%s", defaultPath,expdir,expFname);
datStruct = load(loadExpPath);

bkg_data = permute(squeeze(datStruct.bkg_data(:,:, timeWindow)), [2 1 3]);
measure_data = permute(squeeze(datStruct.measure_data(:,:, timeWindow)), [2 1 3]);

nrows = size(bkg_data,1);
ncols = size(bkg_data,2);
nbinsfinal = ceil(size(bkg_data,3)/expBinWidth);
diff_final = zeros(nrows, ncols, nbinsfinal);

bkg_align = squeeze(bkg_data(alignBkgInds(1), alignBkgInds(2), :));

for r = 1:nrows
    for c = 1:ncols
        
        bkgTpsf_exp = squeeze(bkg_data(r,c,:));
        absTpsf_exp = squeeze(measure_data(r,c,:));

        if (normalizeTransients)
            maxAbs = median(maxk(absTpsf_exp(:), maxk_norm));
            maxBkg = median(maxk(bkgTpsf_exp, maxk_norm));
            absTpsf_exp = absTpsf_exp .* (maxBkg ./ maxAbs);                    % 1. Normalize transients
        end
        
        bkgTpsf_exp = alignTpsf_rising_edge(bkgTpsf_exp, bkg_align, alignInds);
        absTpsf_exp = alignTpsf_rising_edge(absTpsf_exp, bkgTpsf_exp, alignInds);    % 2. Align the transients

        diff = bkgTpsf_exp - absTpsf_exp;
        diff = accumTimeBin(diff, expBinWidth);                                 % 3. Compress data by factor compressFac
        
        diff_final(r, c, :) = diff;
    end
    fprintf("Done with row %d/%d\n", r, nrows);
end

% Plot raw measurements
if (size(diff_final, 3) > 1)
    dispTransImg(diff_final);
end

%% Calculate instrument response function (IRF)
% exp_tpsf = conv(mc_tpsf, irf) -> irf = deconv(exp_tpsf, mc_tpsf)

load(loadBkgPath); % Load background transients

load(loadJpath); % Load Jacobian
vars = fieldnames(Jheaders);
for i = 1:length(vars)
    assignin('base', vars{i}, Jheaders.(vars{i}));
end


% Obtain IRF through deconvolution of experimental bkg data with mc
% simulated bkg data
expBkg = squeeze(bkg_data(bkgDeconvInds(1), bkgDeconvInds(2),:));
expBkg = imgaussfilt(expBkg, imgaussSig);

sizeX = size(expBkg,1)-size(bkgTpsf,1)+1;
fwdMtx = convmtx(bkgTpsf, sizeX);
[f, ft, step] = mat2Handle(fwdMtx, 'mvp');

fprintf("Calculating IRF...");
tic; [irf, ~] = fista(expBkg,f,ft,step,sizeX,irf_retrieval_fistaOpts); irf_retrieval_fistaRuntime = toc;
fprintf("done! Finished in: %.3f seconds\n", irf_retrieval_fistaRuntime);

%% Process simulated Jacobian to match experimental data

lenWindow = length(timeWindow);

% Convolve Jacobian with IRF to match experimental data
J_postproc = ifft(fft(J, lenWindow) .* fft(irf, lenWindow));

% Accumulate rows (corresponding to time bins) of Jacobian
J_postproc = accumTimeBin(J_postproc, expBinWidth);
NBINS = size(J_postproc, 1);

%% Generate full collocated Jacobian

genColocMC_v2;

%% Perform image reconstruction

m_exp = permute(diff_final, [3 1 2]);
m_exp = m_exp(:);

[f_fista, fT_fista, fistaStep] = mat2Handle(J_final, fwd_model);

nVox = VOX_L * VOX_W * VOX_H;
sizeX = [nVox, 1];

lam1Vec = ones(VOX_L, VOX_W, VOX_H);
lam1Vec = reshape(lam1, 1, 1, []) .* lam1Vec;
lam1Vec = lam1Vec(:);

fprintf("Running FISTA...");
tic; [fistaRecon, fistaErr] = fista_layerReg(m_exp,f_fista,fT_fista,fistaStep,sizeX,fistaOpts,lam1Vec); fistaRuntime = toc;
fprintf("done! Finished in: %.3f seconds\n", fistaRuntime);

%% Post-processing (and thresholding)

mu_3D = reshape(fistaRecon, VOX_L, VOX_W, VOX_H);

thresh_recon = zeros(size(mu_3D));
for i = 1:size(thresh_recon, 3)
    recon_i = mu_3D(:,:,i);
    recon_i = recon_i ./ max(recon_i(:));
    thresh_inds = recon_i > threshVec(i);
    recon_i(thresh_inds) = 1;
    recon_i(~thresh_inds) = 0;
    thresh_recon(:,:,i) = recon_i;
end

%% Plot results

layerDx = VOX_W; layerDy = VOX_L; 

% Visualize resulting processed Jacobian
if (visProcessedJ)
    J_reshape = reshape(J_postproc, NBINS, VOX_L, VOX_W, VOX_H);
    figure(); imagesc(squeeze(sum(J_reshape(20:end,:,:,4),1))); % Visualize Jacobian image
    J_plot_trans = squeeze(J_reshape(:,32,32,6)); J_plot_trans = J_plot_trans ./ max(J_plot_trans(:)); % Compare transients with experimental data
    diff_plot_trans = squeeze(diff_final(10,26,:)); diff_plot_trans = diff_plot_trans ./ max(diff_plot_trans(:));
    figure(); plot(J_plot_trans, 'LineWidth', 2); hold on; plot(diff_plot_trans, 'LineWidth', 2);
    legend('Simulated', 'Experimental');
end

if (plotTestIrf)
    figure(); plot(irf);
    convOutLen = 4000;
    convOut = ifft(fft(bkgTpsf, convOutLen) .* fft(irf, convOutLen));
    figure(); hold on; plot(squeeze(bkg_data(16,16,:)), 'lineWidth', 2);  plot(convOut, 'LineWidth', 1);
end

if (plotFistaErr)
    figure(); plot(fistaErr);
end

% 3D Plots
load(truthPath);
truth(truth < 0) = 0; % clamp values to zero for plotting
f_truth = plotRecon3D(truth, 'rgb', sliceImColor, 'cbar', false, 'normalizeIm', false);
view(-158, 38);
colorbar('off');
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'color', 'w');
set(gca, 'FontSize', fontSizeVal);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'visible', 'off');
voxel([SLAB_MINX, SLAB_MINY, SLAB_MINZ],[layerDx, layerDy, layerDz],slab_color,slab_alpha, slab_edge_alpha);
    
mu_3D_gamma = mu_3D .^ recon_gamma;
f_recon = plotRecon3D(mu_3D_gamma, 'rgb', sliceImColor, 'cbar', false, 'normalizeIm', false);
view(-158, 38);
colorbar('off');
set(gcf, 'InvertHardCopy', 'off');
set(gca,'color',[0.96, 0.96, 0.96]);
set(gcf,'color', 'w');
set(gca, 'FontSize', fontSizeVal);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'visible', 'off');
voxel([SLAB_MINX, SLAB_MINY, SLAB_MINZ],[layerDx, layerDy, layerDz],slab_color,slab_alpha, slab_edge_alpha);

f_threshrecon = plotRecon3D(thresh_recon, 'rgb', sliceImColor, 'cbar', false, 'normalizeIm', false);
view(-158, 38);
colorbar('off');
set(gcf, 'InvertHardCopy', 'off');
set(gca,'color',[0.96, 0.96, 0.96]);
set(gcf,'color', 'w');
set(gca, 'FontSize', fontSizeVal);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'visible', 'off');
voxel([SLAB_MINX, SLAB_MINY, SLAB_MINZ],[layerDx, layerDy, layerDz],slab_color,slab_alpha, slab_edge_alpha);

%% Plot Z-slices

% Plot ground truth, reconstructed, and threshold reconstructed images

imH = 2*overall_border + im_border*(length(plotLayers)+1) + VOX_W*length(plotLayers);
imW = 2*overall_border + im_border*(nPlotRows+1) + nPlotRows*VOX_L;
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
    
    truthImTemp = truth(:,:,plotLayers(i)); 
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
    
    reconImTemp = mu_3D(:,:,plotLayers(i));  
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

    threshImTemp = thresh_recon(:,:,plotLayers(i));
    threshImTemp = threshImTemp ./ max(threshImTemp(:));
    rslice = interp1(linspace(0, 1, rgb_axis_len), rax, threshImTemp);
    gslice = interp1(linspace(0, 1, rgb_axis_len), gax, threshImTemp);
    bslice = interp1(linspace(0, 1, rgb_axis_len), bax, threshImTemp);
    threshIm = cat(3, rslice, gslice, bslice);
    slicesIm(threshRow:threshRowEnd, threshCol:threshColEnd, 1:3) = threshIm;

    fprintf("Depth: %f\n", VOX_ORIG(3) + (plotLayers(i)-1)*Jheaders.VOX_ZLEN);
end
fprintf("Full range: %f-%f\n", Jheaders.VOX_ORIGZ, Jheaders.VOX_ORIGZ + (Jheaders.VOX_H * Jheaders.VOX_ZLEN));

f_slices = figure('Position', [155.5000 371 956 552]); imagesc(slicesIm);
axis image
set(gca, 'visible', 'off');

%% Save results

if (shouldSave)
    truth_name = sprintf("%s/%s_%s", savedir, savenameBase, "truth");
    truth_figname = sprintf("%s.fig", truth_name);
    truth_pngname = sprintf("%s.png", truth_name);
    savefig(f_truth, truth_figname);
    export_fig(f_truth, truth_pngname, '-m3', '-transparent', '-png');
    
    recon_name = sprintf("%s/%s_%s", savedir, savenameBase, "recon");
    recon_figname = sprintf("%s.fig", recon_name);
    recon_pngname = sprintf("%s.png", recon_name);
    savefig(f_recon, recon_figname);
    export_fig(f_recon, recon_pngname, '-m3', '-transparent', '-png');
    
    threshrecon_name = sprintf("%s/%s_%s", savedir, savenameBase, "threshrecon");
    threshrecon_figname = sprintf("%s.fig", threshrecon_name);
    threshrecon_pngname = sprintf("%s.png", threshrecon_name);
    savefig(f_threshrecon, threshrecon_figname);
    export_fig(f_threshrecon, threshrecon_pngname, '-m3', '-transparent', '-png');
    
    slices_name = sprintf("%s/%s_%s", savedir, savenameBase, "slices");
    slices_figname = sprintf("%s.fig", slices_name);
    slices_pngname = sprintf("%s.png", slices_name);
    savefig(f_slices, slices_figname);
    export_fig(f_slices, slices_pngname, '-m3', '-transparent', '-png');
end


