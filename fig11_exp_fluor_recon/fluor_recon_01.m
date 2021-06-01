%% Reconstruction parameters

Jdir = "11_28_20_mc_fluor_J";
Jfname = "J";
expDir = "10-Dec-2020_fluor_data";
expFname = "circline2_free-running_coloc_10-Dec-2020_15-53-37";

if (USEGATING)
    binInds = 5400:5800;
else
    binInds = 1:65536;
end

fwdType = "conv";
fistaOpts.lam1 = 4.45e0;
fistaOpts.lam2 = 0.3e-9;
fistaOpts.maxItr = 2000;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;

filts = "E";
nFilts = strlength(filts);
s = 0.0001;

%% Pre-process data

addpath("../lib");

defaultPath = "../dat";
savePath = "../../figs/fig11_exp_fluor_recon";
Jpath = sprintf("%s/%s/%s.mat", defaultPath, Jdir, Jfname);
expPath = sprintf("%s/%s/%s.mat", defaultPath, expDir, expFname);

tic;
load(Jpath);
load(expPath);
loadDatTime = toc;
fprintf("Load data time: %d sec\n", loadDatTime);

tic;
[J_mc, ~] = tempFilts(J, [], Jheaders, filts, 's', s);
tempFiltTime = toc;
fprintf("Temporal filter time: %d sec\n", tempFiltTime);

J_mc = reshape(J_mc, Jheaders.VOX_L, Jheaders.VOX_W);

tic;
[f_fista, fT_fista, fistaStep] = mat2Handle(J_mc, fwdType);
mat2HandleTime = toc;
fprintf("mat2handle time: %d sec\n", mat2HandleTime);

sizeX = [Jheaders.VOX_L, Jheaders.VOX_W];

fluorIm_ungated = sum(squeeze(fluor_data), 3);
fluorIm = sum(squeeze(fluor_data(:,:,:,:,binInds)), 3);

%% Solve inverse problem

fprintf("Running FISTA...");
tic; [fistaRecon, ~] = fista(fluorIm,f_fista,fT_fista,fistaStep,sizeX,fistaOpts); fistaRuntime = toc;
fprintf("done! Finished in: %.3f seconds\n", fistaRuntime);

%% Plot results

fluorIm = fliplr(rot90(fluorIm, 1)); 
fistaRecon = fliplr(rot90(fistaRecon, 1)); 

f1 = figure();
fluorIm_ungated = rot90(fluorIm_ungated,1);
imagesc(fluorIm_ungated);
set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
title("Raw measurements");
axis image;

f2 = figure();
imagesc(fluorIm);
set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
title("Measurements with time-binning");
axis image;

f3 = figure();
imagesc(fistaRecon);
set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
title("Image reconstruction");
axis image;

%% Save results

if ~(exist(savePath, 'dir'))
    mkdir(savePath);
end

if (USEGATING)
    export_fig(f1, sprintf("%s/scene_01_gated_meas.png", savePath), '-m3', '-transparent', '-png');
    export_fig(f2, sprintf("%s/scene_01_gated_time_bins_%d-%d.png", savePath, binInds(1), binInds(end)), '-m3', '-transparent', '-png');
    export_fig(f3, sprintf("%s/scene_01_gated_recon.png", savePath), '-m3', '-transparent', '-png');
else
    export_fig(f1, sprintf("%s/scene_01_ungated_meas.png", savePath), '-m3', '-transparent', '-png');
    export_fig(f3, sprintf("%s/scene_01_ungated_recon.png", savePath), '-m3', '-transparent', '-png');
end

fprintf("Done saving files\n");



