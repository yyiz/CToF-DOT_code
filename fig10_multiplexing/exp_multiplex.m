%% Set up reconstruction parameters

filts = "E";
nFilts = strlength(filts);
s = 0.0001;

plotTacq = [5, 9, 10, 11, 12, 13];
plotRep = 13;

% Automatically set by setting useMultiplex to be true
srcRowInds = 4:7;
srcColInds = 4:7;
detRowInds = 4:7;
detColInds = 4:7;

fwdType = "mvp";
fistaOpts.lam1 = 0.2e-6;
fistaOpts.lam2 = 0.2e-6;
fistaOpts.maxItr = 60;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;

defaultPath = "../dat";

Jfname = sprintf("%s/9_21_20_jacobian_10x10_skull_6p5/J.mat", defaultPath);
bkgFname = sprintf("%s/9_21_20_jacobian_10x10_skull_6p5/tpsf0.mat", defaultPath);
expFname = "intTime_circlepng_full_26-Nov-2020_18-15-20";
expDatFname = sprintf("%s/26-Nov-2020_eink_meas/%s.mat", defaultPath, expFname);
truthFile = sprintf("%s/truth_ims/circle.png", defaultPath);
savedir = "../../figs/fig10_multiplexing_results";
savefname = "EXP_sim_vs_compression_multiplexing";

showTruth = true;

shouldSave = true;

saveVars = ["runFISTA", "runTVAL3", "waveletLvl", "domain", "fistaOpts",...
    "tval3Opts", "admmOpts", "s", "nFilts", "Jfname", "bkgFname", "expDatFname", "truthFile",...
    "savePath", "saveName", "Jheaders", "bkgM_mc", "calibM_mc", "m",...
    "fistaRecon", "fistaReconIm", "fistaErr", "fistaStep", "fistaRuntime"];

%% Pre-process data

addpath("../lib");

% Load (simulated) Jacobian and background TPSF
load(Jfname);
load(bkgFname);

% Load experimental data
load(expDatFname);

selectColoc = true;
useMultiplex = true;
selectSrcDet;

saveName = sprintf("%s/%s.png", savedir, savefname);
saveimdir = sprintf("%s/exp", savedir);

if ~exist(savedir, 'dir')
    mkdir(savedir);
end

if ~exist(saveimdir, 'dir')
    mkdir(saveimdir);
end

%% Set up forward model

[J_mc, bkgM_mc] = tempFilts(J, bkgTpsf, Jheaders, filts, 's', s);

sizeX = size(J_mc,2);
[f_fista, fT_fista, fistaStep] = mat2Handle(J_mc, fwdType,...
    'VOX_L', Jheaders.VOX_L, 'VOX_W', Jheaders.VOX_W);

%% Set up multiplexing

daTruth = rgb2gray(imread(truthFile));
daTruth = double(imresize(daTruth, [Jheaders.VOX_L, Jheaders.VOX_W]));
daTruth = 1 - (daTruth ./ max(daTruth(:)));

psnrVecMultiplex = zeros(length(TacqList), 1);

for intTimeInd = 1:length(TacqList)

bkg_data_temp = bkg_data(:,:,:,:,:,:,intTimeInd);
measure_data_temp = measure_data(:,:,:,:,:,:,intTimeInd);

multiplexExpDatScript;

%%  Calibrate experimental data
nbins = size(bkg_data_temp, 5);
if (nbins > 1)
    bkg_data_temp = sum(bkg_data_temp, 5);
    measure_data_temp = sum(measure_data_temp, 5);
end

calibM_mc = measure_data_temp(:) .* bkgM_mc ./ bkg_data_temp(:);
m = bkgM_mc - calibM_mc;

%% Perform image reconstruction

fprintf("Running FISTA...");
tic; [fistaRecon, fistaErr] = fista(m,f_fista,fT_fista,fistaStep,sizeX,fistaOpts); fistaRuntime = toc;
fprintf("done! Finished in: %.3f seconds\n", fistaRuntime);
  
%% Post-processing and plotting

% Plot FISTA image reconstruction
fistaReconIm = reshape(fistaRecon, [Jheaders.VOX_L, Jheaders.VOX_W, Jheaders.VOX_H]);

if (max(fistaReconIm(:)) > 0)
    fistaReconIm = fistaReconIm ./ max(fistaReconIm(:));
end

psnrVecMultiplex(intTimeInd) = psnr(uint8(255*daTruth), uint8(255*fistaReconIm));

if ismember(intTimeInd, plotTacq)
    multiplexImName = sprintf("%s/multiplex_recon_Tacq=%d.png", saveimdir, TacqList(intTimeInd));
    reconF = figure();
    imagesc(fistaReconIm);
    pbaspect([1 1 1]);
    set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
    colormap(flipud(gray));
    if (shouldSave)
        imwrite(1-fistaReconIm, multiplexImName); % Write image negative
    end
end
end

%% Generate same PSNR plot for unmultiplexed data

% Pre-process data

% Load (simulated) Jacobian and background TPSF
load(Jfname);
load(bkgFname);

% Load experimental data
load(expDatFname);

selectColoc = true;
useMultiplex = false;
selectSrcDet;

% Set up forward model
[J_mc, bkgM_mc] = tempFilts(J, bkgTpsf, Jheaders, filts, 's', s);

sizeX = size(J_mc,2);
[f_fista, fT_fista, fistaStep] = mat2Handle(J_mc, fwdType,...
    'VOX_L', Jheaders.VOX_L, 'VOX_W', Jheaders.VOX_W);

% Set up multiplexing
psnrVec = zeros(length(TacqList), 1);

for intTimeInd = 1:length(TacqList)

    psnrVecTemp = zeros(size(bkg_data, 6), 1);
    
    for repInd = 1:size(bkg_data, 6)
        
        bkg_data_temp = bkg_data(:,:,:,:,:,repInd,intTimeInd);
        measure_data_temp = measure_data(:,:,:,:,:,repInd,intTimeInd);

        %  Calibrate experimental data
        nbins = size(bkg_data_temp, 5);
        if (nbins > 1)
            bkg_data_temp = sum(bkg_data_temp, 5);
            measure_data_temp = sum(measure_data_temp, 5);
        end

        calibM_mc = measure_data_temp(:) .* bkgM_mc ./ bkg_data_temp(:);
        m = bkgM_mc - calibM_mc;

        % Perform image reconstruction

        fprintf("Running FISTA...");
        tic; [fistaRecon, fistaErr] = fista(m,f_fista,fT_fista,fistaStep,sizeX,fistaOpts); fistaRuntime = toc;
        fprintf("done! Finished in: %.3f seconds\n", fistaRuntime);

        % Post-processing and plotting

        fistaReconIm = reshape(fistaRecon, [Jheaders.VOX_L, Jheaders.VOX_W, Jheaders.VOX_H]);

        if (max(fistaReconIm(:)) > 0)
            fistaReconIm = fistaReconIm ./ max(fistaReconIm(:));
        end
        
        psnrVecTemp(repInd) = psnr(uint8(255*daTruth), uint8(255*fistaReconIm));

        if (plotRep == repInd) && ismember(intTimeInd, plotTacq)
            reconName = sprintf("%s/NO_multiplex_recon_Tacq=%d_repInd=%d.png", saveimdir, TacqList(intTimeInd), repInd);
            reconF = figure();
            imagesc(fistaReconIm);
            pbaspect([1 1 1]);
            set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
            colormap(flipud(gray));
            if (shouldSave)
                imwrite(1-fistaReconIm, reconName); % Write image negative
            end
        end
        
    end

    psnrVec(intTimeInd) = mean(psnrVecTemp);

end

%% 

f1 = figure('Position', [680 305 875 673]);
scatter(TacqList, psnrVec, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(TacqList, psnrVecMultiplex, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.5);

interpOrder = 3;
p_nomulti = polyfit(TacqList, psnrVec, interpOrder);
interp_nomulti = polyval(p_nomulti, TacqList);
plot(TacqList, interp_nomulti, 'r', 'LineWidth', 2);

p_multi = polyfit(TacqList, psnrVecMultiplex, interpOrder);
interp_multi = polyval(p_multi, TacqList);
plot(TacqList, interp_multi, 'b', 'LineWidth', 2);

legend("Without multiplexing", "With multiplexing", "Interpolation (without multiplexing)", "Interpolation (with multiplexing)", 'Location', 'NorthEast');

xlabel("Integration time (ms)");
ylabel("PSNR (Image reconstruction quality)");

set(gca, 'XScale', 'log')
set(gca, 'XDir','reverse')
set(gca, 'FontSize', 14);
set(gca, 'FontName', 'Times New Roman');

%% Save results

if (shouldSave)
    if exist(saveName, 'file')
        fprintf("file already exists. Overwrite?\n");
        keyboard;
    end

    saveas(f1, saveName);
    savefig(f1, sprintf("%s/%s", savedir, savefname));
    
    f2 = figure();
    imagesc(daTruth);
    pbaspect([1 1 1]);
    set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
    colormap(flipud(gray));
    
    truthName = sprintf("%s/multiplex_exp_truth.png", saveimdir);
    imwrite(1-daTruth, truthName);
end

