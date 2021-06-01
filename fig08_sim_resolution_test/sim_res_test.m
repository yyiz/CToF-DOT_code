
addpath('../lib');
addpath('../export_fig/');

%% Declare variables

bkg_det_i = 1;
layer_i = 2; % depth of target layer, VOX_ORIG at 5.5
intTime = 5.0; % integration time: seconds
pileupPoint = 5e6; % counts per second

% Fluorescence parameter
useFluor = true; % if true: convolve with fluor lifetime and change noise model
tau = 1e3; % Fluorescence lifetime
od_vec = [0 5]; % Factor determining amount of excitation light reduced by an emission filter;
                

timeBinWidth = 50; % picoseconds
useEndGate = true; % if true: gate from gateStart:end, else use gate
gateStart = 3000;
gate = 1:50;

absCoeff = 0.1;
lineW = 2;
lineLen = 40; 
nreps = 2;

line_scan_row = 64;

% Reconstruction parameters
fistaOpts.lam1 = 0;
fistaOpts.lam2 = 0;
fistaOpts.maxItr = 5e5;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;

% File management parameters
savePath = '../../figs/fig08_sim_resolution_test';
defaultPath = "../dat";
Jdir = "1_19_21_J_128x128_res_test";
bkgDir = "3_24_21_multilayer_128x128vox_bkg";

%% Load in data

if ~exist(savePath, 'dir')
    mkdir(savePath);
end

rng(1);

bkgStruct = load(sprintf("%s/%s/tpsf.mat", defaultPath, bkgDir));
load(sprintf("%s/%s/J.mat", defaultPath, Jdir));
vars = fieldnames(Jheaders);
for i = 1:length(vars)
    assignin('base', vars{i}, Jheaders.(vars{i}));
end

J_cube = reshape(J, NBINS, VOX_L, VOX_W, VOX_H);


J_cube = J_cube(:,:,:,layer_i);
binW = (TIME_MAX - TIME_MIN) / NBINS;
J_mat = reshape(J_cube, NBINS, VOX_L*VOX_W);
[J_postproc, jTimeAx] = convTpsf(J_mat, 'exp', binW, 'lifetime', tau, 'multiCol', true);
NBINS = length(jTimeAx);
Jheaders.NBINS = length(jTimeAx);
J_postproc = reshape(J_postproc, NBINS, VOX_L, VOX_W);

if (useEndGate)
    gate = gateStart:NBINS;
    fprintf("Gate is now: %d-%d\n", gate(1), gate(end));
else
    fprintf("Gate is now: %d-%d\n", gate(1), gate(end));
end

J_gate = J_postproc(gate,:,:);

bkgTpsf = bkgStruct.bkgTpsf(gate, bkg_det_i);
NBINS = length(gate); Jheaders.NBINS = NBINS;

% Compress time bins
J_2D = reshape(J_gate, NBINS, []);
J_2D = accumTimeBin(J_2D, timeBinWidth);
NBINS = size(J_2D, 1); Jheaders.NBINS = NBINS;
J_gate = permute(reshape(J_2D, NBINS, VOX_L, VOX_W), [2 3 1]);
bkg_gate = reshape(accumTimeBin(bkgTpsf, timeBinWidth), 1, 1, []);

timeAx = timeBinWidth .* (0.5:(NBINS-0.5));


% Set up scene
[mu, crow, ccol] = genLinesScene(lineW, VOX_L, VOX_W, absCoeff, 'lineLen', lineLen, 'nreps', nreps);
f_truth = figure(); imagesc(mu);
axis image
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 


export_fig(f_truth, sprintf("%s/ctof-dot_ground_truth.png", savePath), '-transparent', '-m3', '-png');

for r = 1:length(od_vec)
    
reducFac = 10^(-od_vec(r));

%% Generate measurements
fwdType = 'conv';
[f, ft, step] = mat2Handle(J_gate, fwdType, 'nFilts', NBINS);
m = f(mu);

% Apply noise to measurements
maxPhot = intTime * pileupPoint;
bkg_gate = bkg_gate .* reducFac;

m_sig = m + bkg_gate;

m_int = squeeze(sum(m_sig, 3));
maxFac = maxPhot./max(m_int(:));

m_noisy = poissrnd(m_sig .* maxFac)./maxFac;


% Mean filter on m_noisy
timeAx_dim3 = reshape(timeAx, 1, 1, []);
m_noisy = sum(m_noisy .* timeAx_dim3, 3);

% Mean filter on J
J_filt = sum(J_gate .* timeAx_dim3, 3);

[f_fwd, ft_fwd, step_fwd] = mat2Handle(J_filt, fwdType); % Apply mean filter


f_blurred = figure(); imagesc(sum(m_noisy,3));
axis image
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 

sizeX = [VOX_L, VOX_W];
fprintf("Running FISTA...");
tic;
[mu_recon, ~] = fista(m_noisy,f_fwd,ft_fwd,step_fwd,sizeX,fistaOpts);
fista_runtime = toc;
fprintf(sprintf("Done! Runtime: %f sec\n", fista_runtime));
f_recon = figure(); imagesc(mu_recon);
axis image
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
cb = colorbar;
set(findall(f_recon,'-property','FontSize'), 'FontSize', 16);
set(findall(f_recon,'-property','FontName'), 'FontName', 'Times New Roman');

f_line_scan = figure('Position', [680 558 560 243]); plot(mu_recon(line_scan_row, :), 'r', 'LineWidth', 2);
set(findall(f_line_scan,'-property','FontSize'), 'FontSize', 14);
set(findall(f_line_scan,'-property','FontName'), 'FontName', 'Times New Roman');

xtick_vec = get(gca, 'XTick');
xtick_vec = (xtick_vec) * Jheaders.VOX_SIDELEN + Jheaders.VOX_ORIG(1);
xtick_str = sprintfc('%0.1f',xtick_vec);
xticklabels(xtick_str);
%% Calculate MTF

maxVals = zeros(nreps,1);
minVals = zeros(nreps,1);
for i = 1:nreps
    colStart = ccol + (i-1)*2*lineW;
    colEnd = colStart + lineW - 1;
    
    maxArea = mu_recon(crow:(crow+lineLen-1),colStart:colEnd);
    minArea = mu_recon(crow:(crow+lineLen-1),(colStart:colEnd)+lineW);
    
    maxVals(i) = mean(maxArea(:));
    minVals(i) = mean(minArea(:));
end

% Formula based on: https://spie.org/publications/tt52_131_modulation_transfer_function
A_max = mean(maxVals);
A_min = mean(minVals);
M = (A_max-A_min)./(A_max+A_min);

fprintf("Layer: %d, line width: %d, MTF: %f\n", layer_i, lineW, M);

%% Save results


savename = sprintf("ctof-dot_recon_od=%d_gate=%d-%d", od_vec(r), gate(1), gate(end));
strrep(savename, '.', 'p');
export_fig(f_blurred, sprintf("%s/%s_meas.png", savePath, savename), '-transparent', '-m3', '-png');
export_fig(f_recon, sprintf("%s/%s_recon.png", savePath, savename), '-transparent', '-m3', '-png');
export_fig(f_line_scan, sprintf("%s/%s_line_scan.png", savePath, savename), '-transparent', '-m3', '-png');
save(sprintf("%s/%s_params", savePath, savename), "fistaOpts",...
    "layer_i", "lineW", "intTime", "pileupPoint", "tau", "useEndGate",...
    "gateStart", "gate", "absCoeff", "lineLen", "nreps", "VOX_L", "VOX_W",...
    "Jheaders", "bkg_gate", "J_gate", "mu", "crow", "ccol", "m", "m_noisy", "mu_recon",...
    "A_max", "A_min", "M", "fista_runtime");

end

