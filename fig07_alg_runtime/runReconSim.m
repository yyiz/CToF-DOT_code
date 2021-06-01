% clear; close all; clc;

rng(1);

% Assumes Jacobian of correct dimensions has already been generated
truthFname = "R2";

bw_thresh = 60; % When converting png to binary image, this is pixel thresholding (0-255)

delUa = 0.1;

% Noise model
intTime = 1e0; % seconds
maxPhot = intTime * 5e6;

% Time-gating window

if (intSig)
    timebinInds = 1:10; % Select early arriving photons
else
    timebinInds = selectedBins;
end

% Reconstruction parameters
fistaOpts.lam1 = 0;
fistaOpts.lam2 = 0;
fistaOpts.maxItr = 100;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;

%% Load in Jacobian

vars = fieldnames(Jheaders);
for i = 1:length(vars)
    assignin('base', vars{i}, Jheaders.(vars{i}));
end

if (startsWith(fwd_type, 'mvp'))
    % If MVP, do nothing
elseif (startsWith(fwd_type, 'conv'))
    J = reshape(J, [NBINS, VOX_W, VOX_L]);
    J = permute(J, [2, 3, 1]);   
else
    fprintf("Invalid forward model type\n");
end

%% Set up scene (Generate mu vector)

truthIm = sprintf("truth_ims/%s.png", truthFname);
mu = imread(sprintf("../dat/%s", truthIm));
mu = imresize(mu, [VOX_L, VOX_W]);
if (size(mu,3) > 1)
    mu = rgb2gray(mu);
end
mu = double(mu);
absInds = mu < bw_thresh;
mu(absInds) = delUa;
mu(~absInds) = 0;    

if (startsWith(fwd_type, 'mvp'))
    mu_input = mu(:);
elseif (startsWith(fwd_type, 'conv'))
    mu_input = repmat(mu, [1 1 NBINS]);
else
    fprintf("Invalid forward model type\n");
end

%% Generate differential measurements

addpath("../lib/");
[f, ~, ~] = mat2Handle(J, fwd_type,...
                       'VOX_L', VOX_L, 'VOX_W', VOX_W);

m = f(mu_input);
m(m < 0) = 0; % Remove negligible negative values, otherwise poissrnd will result in nan

%% Apply time-gating

if (startsWith(fwd_type, 'mvp'))
    NVOX = Jheaders.VOX_L * Jheaders.VOX_W * Jheaders.VOX_H;
    J_binned = reshape(J, NBINS, [], NVOX);
    m_binned = reshape(m, NBINS, []);
    
    if (intSig)
        J_binned = sum(J_binned(timebinInds, :, :), 1);
        m_binned = sum(m_binned(timebinInds, :, :), 1);
        bkg_binned = sum(bkgTpsf(timebinInds, :), 1);
        NBINS = 1;        
    else
        J_binned = J_binned(timebinInds, :, :);
        m_binned = m_binned(timebinInds, :, :);
        bkg_binned = bkgTpsf(timebinInds, :);
        NBINS = length(timebinInds);
    end
    
    J_binned = reshape(J_binned, [], NVOX);
    m_binned = m_binned(:);
    bkg_binned = bkg_binned(:);
    
elseif (startsWith(fwd_type, 'conv'))
    J_binned = J(:,:,timebinInds);
    m_binned = m(:,:,timebinInds);

    bkg_binned = permute(bkgTpsf(timebinInds), [2 3 1]);
    
    NBINS = length(timebinInds);
else
    fprintf("Invalid forward model type\n");
end

%% Apply noise
    
if (startsWith(fwd_type, 'mvp'))
    maxMeas = sum(reshape(bkg_binned, NBINS, []), 1);
    maxMeas = max(maxMeas(:));
else
    maxMeas = sum(bkg_binned(:)); % There should only be 1 TPSF in conv mode
end   

normFac = maxPhot / maxMeas;

absMeas = abs(bkg_binned - m_binned);

bkg_binned = bkg_binned .* normFac;
absMeas = absMeas .* normFac;
m_binned = poissrnd(bkg_binned) - poissrnd(absMeas);
m_binned = m_binned ./ normFac;


%% Run reconstruction

[f, fT, step] = mat2Handle(J_binned, fwd_type,...
                            'VOX_L', VOX_L, 'VOX_W', VOX_W);

if (startsWith(fwd_type, 'mvp'))
    sizeX = NVOX;
else
    sizeX = size(mu);
end
                        
fprintf("Running FISTA...");
tic; [reconFISTA, ~] = fista(m_binned, f, fT, step, sizeX, fistaOpts); fistaRuntime = toc;
fprintf("done! Finished in: %.3f seconds\n", fistaRuntime);

%% Plot results

if (startsWith(fwd_type, 'mvp'))
    reconFISTA = reshape(reconFISTA, VOX_L, VOX_W, VOX_H);
else
    % Do not need to reshape output image if running deconvolution
end

f1 = figure('Position', [34.5000 580 1035 348.5000]);

subplot(1,2,1); imagesc(mu); title("Ground truth");
set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
subplot(1,2,2);
imagesc(reconFISTA);
set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
title("Reconstruction: time-gated");

set(findall(f1,'-property','FontSize'), 'FontSize', 14);
set(findall(f1,'-property','FontName'), 'FontName', 'Times New Roman');

