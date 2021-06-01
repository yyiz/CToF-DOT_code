clear; close all; clc;

defaultPath = "../dat";
mcDatDir = "11_8_20_mc_optical_params";
mcDatFname = "mus_tpsf_mat";

plotInds = 4800:8000; %Focus on TPSF region
tpsfInd = 3; % Which MC TPSF to compare to, the scattering coefficient is obtained through `mus_vec(tpsfInd)`
measBinSize = 1; % Binwidth for experimental data, units: ps
shiftAmt = 120; % Amount to shift measured TPSF when aligning to reference

%% Compare measured TPSF with MC

load(sprintf("%s/%s/%s", defaultPath, mcDatDir, mcDatFname));

% Make sure MC simulation and experimental data are aligned to the same time axes using interpolation
DTIME = (mHeaders.TIME_MAX - mHeaders.TIME_MIN)/mHeaders.NBINS;
musTpsf = tpsf_mat(tpsfInd,:);
sampTimeAx = timeAx .* (measBinSize / DTIME);
sampMusTpsf = interp1(timeAx, musTpsf, sampTimeAx);
sampMusTpsf = sampMusTpsf(~isnan(sampMusTpsf));

measDir = "11_10_20_mc_optical_params";
measFname = "11_10_20_retrieve_optical_properties_3cmX3cmX2cm";

load(sprintf("%s/%s/%s", defaultPath, measDir, measFname));

irfTpsf = double(irf');
measTpsf = double(meas');

% Convolve MC result with IRF
convTpsf = conv(irfTpsf, sampMusTpsf);

% Normalize TPSF amplitudes (account for differences due to laser power,
% etc., only want to test effects of scattering)
convTpsf = convTpsf ./ max(convTpsf(:));
measTpsf = measTpsf ./ max(measTpsf(:));

% Shift measTpsf to align with reference TPSF
measTpsf = circshift(measTpsf, shiftAmt);

f1 = figure();
plot(convTpsf(plotInds), 'LineWidth', 2);
hold on;
plot(measTpsf(plotInds), 'LineWidth', 2);
legend('MC Output', 'Measured');
set(gca, 'FontSize', 18);
set(gca, 'FontName', 'Times New Roman');
fprintf("Testing scattering properties for mu_s=%d\n", mus_vec(tpsfInd));


savefname = sprintf("match_tpsf_mus=%d", mus_vec(tpsfInd));
savedir = "../../figs/fig09_exp_resolution_test";
saveName = sprintf("%s/%s.png", savedir, savefname);

if ~exist(savedir, 'dir')
    mkdir(savedir);
end

if exist(saveName, 'file')
    fprintf("file already exists. Overwrite?\n");
    keyboard;
end
saveas(f1, saveName);