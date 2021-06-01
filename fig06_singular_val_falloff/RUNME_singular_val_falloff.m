clear; close all; clc;

%%

fdir = "12_10_20_Jacobian_5x5_singular_val";
fname = "J_srcdet_5x5";
savedir = "../../figs/fig06_singular_val_falloff";
savename = sprintf("%s_svd", fname);

fullLoadName = sprintf("../dat/%s/%s", fdir, fname);
tic;
load(fullLoadName);
toc;

nvox = Jheaders.VOX_L * Jheaders.VOX_W * Jheaders.VOX_H;

%% Generate different Jacobians

addpath("../lib");
addpath("../export_fig");

timeGateInds = 3:Jheaders.NBINS;
J_time_dim1 = reshape(J, Jheaders.NBINS, []);

J_int = sum(J_time_dim1, 1);
J_int = reshape(J_int, [], nvox);

tic; S_allbins = svd(J); t1 = toc;
fprintf("SVD(J) time: %d\n", t1);
S_int = svd(J_int);

% Normalize singular value plots
S_allbins = S_allbins ./ max(S_allbins(:));
S_int = S_int ./ max(S_int(:));

% Load in collocated 25 x 25 
load("../dat/12_12_20_coloc_Jacobian_5x5_singular_val/singular_vals_J_coloc.mat");

% Condition numbers
K_allbins = max(S_allbins) ./ min(S_allbins);
K_int = max(S_int) ./ min(S_int);
K_coloc_25x25 = max(S_coloc_25x25) ./ min(S_coloc_25x25);
%% Plot results

f1 = figure('Position', [251 140.5000 1326 817]); 

semilogy(S_int, 'LineWidth', 5);
hold on;
semilogy(S_allbins, 'LineWidth', 5);
semilogy(S_coloc_25x25, 'LineWidth', 5);

ylim([1e-3, 1e0]);
set(gca, 'FontSize', 30);
set(gca, 'FontName', 'Times New Roman');

noiseThresh = 1e-2;

max_ind_int = find(S_int >= noiseThresh, 1, 'last');
max_ind_allbins = find(S_allbins >= noiseThresh, 1, 'last');
max_ind_coloc_25x25 = find(S_coloc_25x25 >= noiseThresh, 1, 'last');

fprintf("Index of minimum singular value greater than 1e-3:\n\n");
fprintf("Traditional DOT: %d\n", max_ind_int);
fprintf("All time bins: %d\n", max_ind_allbins);
fprintf("Confocal: %d\n", max_ind_coloc_25x25);
%% Save results

savepath = sprintf("%s/%s", savedir, savename); 
matname = sprintf("%s.mat", savepath);

if (exist(matname, 'file'))
    fprintf("File exists, save anyway?\n");
    keyboard;
end
if ~(exist(savedir, 'dir'))
    mkdir(savedir);
end

save(matname, "S_int", "S_coloc_25x25",...
    "S_allbins", "timeGateInds", "Jheaders", "fullLoadName");
savefig(f1, savepath);
export_fig(f1, sprintf("%s.png", savepath), '-m3', '-transparent', '-png');
fprintf("Done saving files\n");
