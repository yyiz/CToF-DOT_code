
Jdir = "3_21_21_sim_3D_Jacobian";
Jfname = "J_coloc_nSD=32x32_voxSize=32x32x6";

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
savedir = sprintf("%s/%s", savedir, savename);

plotScene = true;

%% Generate and visualize scene
addpath("../lib");
addpath("../export_fig");

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
plotMu = mu;

layerDx = VOX_W; layerDy = VOX_L; 

% Plot scene
if (plotScene)
    ftruth = plotRecon3D(plotMu, 'rgb', sliceImColor, 'cbar', false, 'normalizeIm', false);
    view([227 35]);
    set(ftruth, 'InvertHardCopy', 'off');
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); set(gca,'ZTickLabel',[]);
    xlabel(''); ylabel(''); zlabel('');

    cb = colorbar;
    maxPix = max(plotMu(:));
    cb.Ticks = linspace(0,maxPix,numTicks_truth);
    tlabels = [];
    for i = 1:numTicks_truth
        tlabels = [tlabels sprintf(cbarFmt, cb.Ticks(i))];
    end
    cb.TickLabels = tlabels;
    pos =  cb.Position; %gets the positon and size of the color bar
    set(cb,'Position',[pos(1)+dx pos(2)+dy pos(3)+dw pos(4)+dh])% To change size
    set(gca,'visible', 'off');
    set(gcf,'color', 'w');
    
    set(gca, 'FontSize', fontSizeVal);
    set(gca, 'FontName', 'Times New Roman');
    
    voxel([SLAB_MINX, SLAB_MINY, SLAB_MINZ],[layerDx, layerDy, layerDz],slab_color,slab_alpha, slab_edge_alpha);
end

mu = mu(:);

if (shouldSave)
    truthFigPath = sprintf("%s/truth_%s.fig", savedir, savename);
    truthPngPath = sprintf("%s/truth_%s.png", savedir, savename);
    if ~(exist(savedir, 'dir'))
        mkdir(savedir);
    end
    
    savefig(ftruth, truthFigPath);
    export_fig(ftruth, truthPngPath, '-m3', '-transparent', '-png');
end

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

% Accumulate time bins
J = reshape(J, NBINS, nSrcDet*nVox);

J = J(selectBins,:);
bkgTpsf_binned = bkgTpsf(selectBins);
NBINS = length(selectBins);
% Form final Jacobian
J = reshape(J, NBINS*nSrcDet, nVox);

m = J*mu;

% Apply poisson noise
if (useFluor)
    maxIntM = max(sum(reshape(m, Jheaders.NBINS, []), 1));
    scaleFac = maxPhots / maxIntM;
    m = poissrnd(m .* scaleFac);
    m = m ./ scaleFac;
else
    m_reshape = reshape(m, NBINS, nSrcDet);
    intBkgSig = sum(bkgTpsf_binned(:));
    scaleFac = maxPhots ./ intBkgSig;
    
    bkgMeas = repmat(bkgTpsf_binned, [1, nSrcDet]);
    absMeas = abs(bkgMeas - m_reshape); % Handle numerical imprecision with abs
    
    bkgMeas = bkgMeas .* scaleFac;
    absMeas = absMeas .* scaleFac;
    m = poissrnd(bkgMeas) - poissrnd(absMeas);
    
    m = m ./ scaleFac;
    m = m(:);
end

%% Solve inverse problem

[f, fT, step] = mat2Handle(J, fwd_type);

nVox = VOX_L * VOX_W * VOX_H;
sizeX = [nVox, 1];

lam1Vec = ones(VOX_L, VOX_W, VOX_H);
lam1Vec = reshape(lam1, 1, 1, []) .* lam1Vec;
lam1Vec = lam1Vec(:);

fprintf("Running FISTA...");
tic; [fistaRecon, ~] = fista_layerReg(m,f,fT,step,sizeX,fistaOpts,lam1Vec); fistaRuntime = toc;
fprintf("done! Finished in: %.3f seconds\n", fistaRuntime);

fistaRecon = reshape(fistaRecon, [Jheaders.VOX_W, Jheaders.VOX_L, Jheaders.VOX_H]);

%% Plot results

VOX_ORIG = [Jheaders.VOX_ORIGX, Jheaders.VOX_ORIGY, Jheaders.VOX_ORIGZ];
VOXDIM = [Jheaders.VOX_SIDELEN, Jheaders.VOX_SIDELEN, Jheaders.VOX_SIDELEN];
fistaRecon_hdr = fistaRecon .^ gamma;
frecon = plotRecon3D(fistaRecon_hdr, 'rgb', sliceImColor, 'cbar', false, 'normalizeIm', false);
view([227 35]);
set(frecon, 'InvertHardCopy', 'off');
xlabel(''); ylabel(''); zlabel('');

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
set(gcf,'color', 'w');
set(gca, 'FontSize', fontSizeVal);
set(gca, 'FontName', 'Times New Roman');
set(gca,'visible', 'off');
voxel([SLAB_MINX, SLAB_MINY, SLAB_MINZ],[layerDx, layerDy, layerDz],slab_color,slab_alpha, slab_edge_alpha);


reconFigPath = sprintf("%s/recon_%s.fig", savedir, savename);
reconPngPath = sprintf("%s/recon_%s.png", savedir, savename);

if (shouldSave)
    savefig(frecon, reconFigPath);
    export_fig(frecon, reconPngPath, '-m3', '-transparent', '-png');
end
%% Display threshold image

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

fthresh = plotRecon3D(thresh, 'rgb', sliceImColor);
view([227 35]);
set(fthresh, 'InvertHardCopy', 'off');
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); set(gca,'ZTickLabel',[]);
set(gcf,'color', 'w');
set(gca, 'FontSize', fontSizeVal);
set(gca, 'FontName', 'Times New Roman');
xlabel(''); ylabel(''); zlabel('');
colorbar('off');
set(gca,'visible', 'off');
voxel([SLAB_MINX, SLAB_MINY, SLAB_MINZ],[layerDx, layerDy, layerDz],slab_color,slab_alpha, slab_edge_alpha);


threshFigPath = sprintf("%s/thresh_%s.fig", savedir, savename);
threshPngPath = sprintf("%s/thresh_%s.png", savedir, savename);
if (shouldSave)
    savefig(fthresh, threshFigPath);
    export_fig(fthresh, threshPngPath, '-m3', '-transparent', '-png');
end

%% Plot Z-slices

% Plot ground truth, reconstructed, and threshold reconstructed images
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
    
    truthImTemp = plotMu(:,:,plotLayers(i)); 
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

slicesFigPath = sprintf("%s/slices_%s.fig", savedir, savename);
slicesPngPath = sprintf("%s/slices_%s.png", savedir, savename);

f_slices = figure('Position', [155.5 107.5 1184 815.5]); imagesc(slicesIm);
axis image
set(gca, 'visible', 'off');

if (shouldSave)
    savefig(f_slices, slicesFigPath);
    export_fig(f_slices, slicesPngPath, '-m3', '-transparent', '-png');
end






