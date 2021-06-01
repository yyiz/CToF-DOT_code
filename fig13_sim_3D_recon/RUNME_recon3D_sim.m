clear; close all; clc;

%% Artifical vasculature scene

truthImName = "truth_vasculature_32x32x6";

gamma = 0.6;
numTicks_recon = 6;
numTicks_truth = 2;

absMua = 0.01 .* [1 1 1 1 1 1];
intTime = 1e1;
maxPhots = 5e6 * intTime;

selectBins = 1:10;

useFluor = false;
tau = 1000;

fwd_type = 'mvp';
fistaOpts.lam2 = 0;
fistaOpts.maxItr = 175;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;
lam1 = [2e-7, 2e-8, 0, 0, 0, 0];

threshVec = [0.2, 0.02, 0.4, 0.4, 0.7, 0.7];

recon3D_sim;

%% Circle-lines scene

truthImName = "truth_circ_32x32x6";

gamma = 0.6;
numTicks_recon = 6;
numTicks_truth = 7;

absMua = 0.01 .* [1 2 3 4 5 6];
intTime = 1e1;
maxPhots = 5e6 * intTime;

selectBins = 1:10;

useFluor = false;
tau = 1000;

fwd_type = 'mvp';
fistaOpts.lam2 = 0;
fistaOpts.maxItr = 75;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;
lam1 = [7.2e-6, 19e-7, 95e-8, 60e-8, 15e-8, 80e-9];

threshVec = [0.09, 0.2, 0.1, 0.1, 0.2, 0.2];

recon3D_sim;

%% Deconvolution circle-line scene

truthImName = "truth_circ_32x32x6";
deconv3D;

