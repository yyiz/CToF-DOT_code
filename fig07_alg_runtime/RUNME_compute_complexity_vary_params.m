clear; close all; clc;

addpath("../lib");

%% Vary voxel size

previewScene = false;
basedir = "../../figs/fig07_alg_runtimes";
savedir = sprintf("%s/vox_size", basedir);

if ~exist(savedir, 'dir')
    mkdir(savedir);
end

voxDims = [16, 20, 32, 40, 64];
myLegend = ["DOT", "ToF-DOT", "Ours"];
myLegendInterp = ["Interp (Traditional)", "Interp (ToF-DOT)", "Interp (Ours)"];
runtimeArr = zeros(length(voxDims), length(myLegend));
psnrArr = zeros(length(voxDims), length(myLegend));

% If using timebinning, select the bin to use
selectedBins = 5:10;

% Source path parameters (mm)
SRC_ORIG = [-26.5; -26.5; 0.0];
SRC_L = 8;
SRC_W = 8;
SRC_SEP = 7.5;
% Detector path parameters (mm)
DET_ORIG = [-26.5; -26.5; 0.0];
DET_L = 8;
DET_W = 8;
DET_SEP = 7.5;
DET_LEN = 1;

% Forward model to employ
fwd_type = 'mvp';

for r = 1:length(myLegend)
    
if (r == 1)
    intSig = true;
    coloc = false;
elseif (r == 2)
    intSig = false;
    coloc = false;
else
    intSig = false;
    coloc = true;
end
    
for k = 1:length(voxDims)

    voxDim = voxDims(k);
    
    VOX_ORIG = [-32; -32; 6.5];
    VOX_W = voxDim;
    VOX_L = voxDim;
    VOX_H = 1;
    VOXDIM = [64/voxDim; 64/voxDim; 1.0];

    % Compute Jacobian
    compute_complexity_vary_params;

    % Perform reconstruction
    runReconSim
    
    runtimeArr(k, r) = fistaRuntime;
    
    truthIm = uint8(255 .* (mu ./ max(mu(:))));
    fistaIm = uint8(255 .* (reconFISTA ./ max(reconFISTA(:))));
    psnrArr(k, r) = psnr(truthIm, fistaIm);
    
    fname = sprintf("voxdim_%d_alg_%s.png", voxDim, myLegend(r));
    imwrite(255-fistaIm, sprintf("%s/%s", savedir, fname));
end
end

interpOrder = 2;
f1_vox = figure('Position', [590.5000 489 435 406.5000]);
hold on;
h = [];
for i = 1:length(myLegend)
    if (i == 1)
        dotColor = 'g';
    elseif (i == 2)
        dotColor = 'b';
    else
        dotColor = 'r';
    end
    plotI = runtimeArr(:,i);
    scatter(voxDims, plotI, 10, dotColor, 'filled', 'MarkerFaceAlpha', 0.5);
    
    p = polyfit(voxDims', plotI, interpOrder);
    interpP = polyval(p, voxDims);
    h = [h plot(voxDims, interpP, dotColor, 'LineWidth', 2,...
        'DisplayName', myLegend(i))];
end

set(gca, 'FontSize', 16);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'YScale', 'log')

xtick_vec = get(gca, 'XTick');
xtick_str = sprintfc('%d',xtick_vec.^2);
xticklabels(xtick_str);


f2_vox = figure(); plot(voxDims, psnrArr, 'LineWidth', 2);
legend(myLegend, 'Location', 'NorthEast');

if (exist(sprintf("%s/runtimes_voxsize.png", basedir), 'file'))
    fprintf("File exists. Overwrite?\n");
    keyboard;
end

saveas(f1_vox, sprintf("%s/runtimes_voxsize.png", basedir));
savefig(f1_vox, sprintf("%s/runtimes_voxsize.fig", basedir));
saveas(f2_vox, sprintf("%s/psnr_voxsize.png", basedir));
savefig(f2_vox, sprintf("%s/psnr_voxsize.fig", basedir));
save(sprintf("%s/runtimes_voxdims", basedir), "voxDims", "psnrArr", "runtimeArr");


%% Vary source-detector array size

previewScene = false;
basedir = "../../figs/fig07_alg_runtimes";
savedir = sprintf("%s/nsrc_det", basedir);

if ~exist(savedir, 'dir')
    mkdir(savedir);
end

nSrcDet = 4:10;

myLegend = ["DOT", "ToF-DOT", "Ours"];
myLegendInterp = ["Interp (Traditional)", "Interp (ToF-DOT)", "Interp (Ours)"];

runtimeArr = zeros(length(nSrcDet), length(myLegend));
psnrArr = zeros(length(nSrcDet), length(myLegend));

% Select time bins to use
selectedBins = 5:10;

% Forward model to employ
fwd_type = 'mvp';

for r = 1:length(myLegend)
    
if (r == 1)
    intSig = true;
    coloc = false;
elseif (r == 2)
    intSig = false;
    coloc = false;
else
    intSig = false;
    coloc = true;
end
    
for k = 1:length(nSrcDet)
    % Jacobian voxel parameters
    VOX_ORIG = [-32; -32; 6.5];
    VOX_W = 32;
    VOX_L = 32;
    VOX_H = 1;
    VOXDIM = [2.0; 2.0; 1.0];

    SRC_ORIG = [-26.5; -26.5; 0.0];
    DET_ORIG = [-26.5; -26.5; 0.0];
    DET_LEN = 1;
    
    
    nSD = nSrcDet(k);
    
    % Source path parameters (mm)
    SRC_L = nSD;
    SRC_W = nSD;
    SRC_SEP = 7.5 * (8/nSD);
    % Detector path parameters (mm)
    DET_L = nSD;
    DET_W = nSD;
    DET_SEP = 7.5 * (8/nSD);
    
    % Compute Jacobian
    compute_complexity_vary_params;

    % Perform reconstruction
    runReconSim
    
    runtimeArr(k, r) = fistaRuntime;
    
    truthIm = uint8(255 .* (mu ./ max(mu(:))));
    fistaIm = uint8(255 .* (reconFISTA ./ max(reconFISTA(:))));
    psnrArr(k, r) = psnr(truthIm, fistaIm);

    fname = sprintf("nsrcdet_%d_alg_%s.png", nSD, myLegend(r));
    imwrite(255-fistaIm, sprintf("%s/%s", savedir, fname));
end
end

interpOrder = 2;
f1_srcdet = figure('Position', [590.5000 489 435 406.5000]);
hold on;
h = [];
for i = 1:length(myLegend)
    if (i == 1)
        dotColor = 'g';
    elseif (i == 2)
        dotColor = 'b';
    else
        dotColor = 'r';
    end
    plotI = runtimeArr(:,i);
    scatter(nSrcDet, plotI, 10, dotColor, 'filled', 'MarkerFaceAlpha', 0.5);
    
    p = polyfit(nSrcDet', plotI, interpOrder);
    interpP = polyval(p, nSrcDet);
    h = [h plot(nSrcDet, interpP, dotColor, 'LineWidth', 2,...
        'DisplayName', myLegend(i))];
end

set(gca, 'FontSize', 16);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'YScale', 'log')

xtick_vec = get(gca, 'XTick');
xtick_vec = reshape([xtick_vec; xtick_vec], 1,[]);
xtick_str = sprintfc('%dx%d',xtick_vec);
xticklabels(xtick_str);

f2_srcdet = figure(); plot(nSrcDet, psnrArr, 'LineWidth', 2);
legend(myLegend, 'Location', 'NorthEast');

if (exist(sprintf("%s/runtimes_nsrcdet.png", basedir), 'file'))
    fprintf("File exists. Overwrite?\n");
    keyboard;
end

saveas(f1_srcdet, sprintf("%s/runtimes_nsrcdet.png", basedir));
savefig(f1_srcdet, sprintf("%s/runtimes_nsrcdet.fig", basedir));
saveas(f2_srcdet, sprintf("%s/psnr_nsrcdet.png", basedir));
savefig(f2_srcdet, sprintf("%s/psnr_nsrcdet.fig", basedir));
save(sprintf("%s/runtimes_nsrcdet", basedir), "nSrcDet", "psnrArr", "runtimeArr");

