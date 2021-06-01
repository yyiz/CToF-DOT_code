clear; close all; clc;

%% Reconstruction of 01 target scene without gating

USEGATING = false;
fluor_recon_01;

%% Reconstruction of 01 target scene with gating

USEGATING = true;
fluor_recon_01;

%% Reconstruction of 2 lines scene

fluor_recon_2lines