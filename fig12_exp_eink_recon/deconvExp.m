%% Parameters

measType = 'conv';
s = 0.001;
filts = "M";
nFilts = strlength(filts);
filtWeights = [1 1 1];

%% Pre-process data

measDir = "02-Dec-2020_eink_meas";
measFilename = sprintf("%s/%s", measDir, measFile);
jacobianFilename = "11_11_20_deconv_exp_psf/J";

savepath = "../../figs";
saveDir = "fig12_exp_eink_recon";

binSz = 1;
timeGate = 1:65536;

isFluorDat = false;

defaultPath = "../dat";
savePath = sprintf("%s/%s", savepath, saveDir);
saveName = sprintf("recon_%s", measFile);

addpath("../lib/");

load(sprintf("%s/%s", defaultPath, jacobianFilename));
VOX_L = Jheaders.VOX_L; VOX_W = Jheaders.VOX_W;

truthIm = imread(sprintf("%s/%s", defaultPath, truthImFile));
loadMeas = load(sprintf("%s/%s", defaultPath, measFilename));
if (isFluorDat)
    nMeasBins = size(loadMeas.fluor_data, 5);
    m = loadMeas.fluor_data(:,:,:,:,1:nMeasBins);
else
    loadMeas.bkg_data = loadMeas.bkg_data(:,:,:,:,timeGate);
    loadMeas.measure_data = loadMeas.measure_data(:,:,:,:,timeGate);
    nMeasBins = size(loadMeas.bkg_data, 5);
    bkg = loadMeas.bkg_data(:,:,:,:,1:nMeasBins);
    abs = loadMeas.measure_data(:,:,:,:,1:nMeasBins);
    m = bkg - abs;
end
m = squeeze(permute(m, [5 1 2 3 4]));
m = reshape(m, nMeasBins, []);

Jheaders.nMeasBins = nMeasBins;
Jheaders.measTimeAx = transpose(binSz:binSz:(binSz*nMeasBins));
[J, m] = tempFilts(J, m, Jheaders, filts, 's', s, 'filtWeights', filtWeights, 'fwdModelType', measType);

if (strcmp(measType, "conv"))
    sizeX = [VOX_L, VOX_W];
    J = reshape(J, VOX_W, VOX_L, nFilts);
    m = reshape(m, VOX_W, VOX_L, nFilts);
else
    sizeX = [VOX_L*VOX_W, 1];
end

J = J ./ max(J(:));
m = m ./ max(m(:));

[f, fT, step] = mat2Handle(J, measType);

%%

fprintf("Running FISTA...");
tic; [reconFista, fistaErr] = fista(m,f,fT,step,sizeX,fistaOpts); fistaRuntime = toc;
fprintf("done! Runtime: %0.3f sec\n", fistaRuntime);

reconFista = reshape(reconFista, Jheaders.VOX_W, Jheaders.VOX_L);
reconFista = rot90(reconFista, 1);
m = rot90(m, 1);

f1 = figure();
imagesc(truthIm);
axis image;
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

f2 = figure();
imagesc(m);
axis image;
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

f3 = figure();
imagesc(reconFista);
axis image;
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

%% Save results

if ~(exist(savePath, 'dir'))
    mkdir(savePath);
end

export_fig(f1, sprintf("%s/deconv_%s_truth.png", savePath, measFile), '-m3', '-transparent', '-png');
export_fig(f2, sprintf("%s/deconv_%s_meas.png", savePath, measFile), '-m3', '-transparent', '-png');
export_fig(f3, sprintf("%s/deconv_%s_recon.png", savePath, measFile), '-m3', '-transparent', '-png');




