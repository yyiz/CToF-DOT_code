%% Set up reconstruction parameters

runFISTA = true;

fistaOpts.lam1 = 0;
fistaOpts.lam2 = 0;
fistaOpts.maxItr = 250;
fistaOpts.tol = 1e-12;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;

defaultPath = "../dat";
measFilename = "fluor_16-Mar-2020/fluor_coloc_16-Mar-2020_20-29-35.mat";
jacobianFilename = "6_4_20_reconFluorExp/jacobian0.mat";

saveVars = ["fistaOpts", "shouldSave", "savePath", "saveName",...
    "m", "psf"];

shouldSave = true;
savePath = "../../figs/fig11_exp_fluor_recon";
saveName = "recon_2lines";
    
%% Pre-process data
    
addpath("../lib/");

loadMeas = load(sprintf("%s/%s", defaultPath, measFilename));
load(sprintf("%s/%s", defaultPath, jacobianFilename));

m = sum(squeeze(loadMeas.fluor_data), 3);
m = flipud(transpose(m));
sizeX = [params.VOX_L, params.VOX_W];

psf = sum(J, 1);
psf = reshape(psf, params.VOX_L, params.VOX_W);
psf = psf';

[f, fT, step] = mat2Handle(psf, "conv");

%% Perform image reconstruction

fprintf("Running FISTA...");
tic; [reconFista, fistaErr] = fista(m,f,fT,step,sizeX,fistaOpts); fistaRuntime = toc;
fprintf("done! Runtime: %0.3f sec\n", fistaRuntime);

%% Plot results

f1 = figure();
imagesc(m);
pbaspect([1 1 1]);
title("Measurements");
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

f2 = figure();
imagesc(reconFista);
pbaspect([1 1 1]);
title("FISTA");
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

%% Save results

if ~exist(savePath, 'dir')
    mkdir(savePath);
end

saveVars = [saveVars, "fistaErr", "reconFista", "fistaRuntime"];

matName = sprintf("%s/%s.mat", savePath, saveName);
saveVars = cellstr(saveVars);
save(matName, "-v7.3", saveVars{:});

export_fig(f1, sprintf("%s/2lines_meas.png", savePath), '-m3', '-transparent', '-png');
export_fig(f2, sprintf("%s/2lines_recon.png", savePath), '-m3', '-transparent', '-png');

