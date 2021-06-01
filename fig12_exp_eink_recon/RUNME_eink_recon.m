clear; close all; clc;

shouldSave = true;

% Reconstruct circline2
measFile = "circline2png_coloc_02-Dec-2020_21-27-34";
truthImFile = "truth_ims/circline2.png";
fistaOpts.lam1 = 5e0;
fistaOpts.lam2 = 0;
fistaOpts.maxItr = 50;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;

deconvExp;

%% 

clear;

shouldSave = true;

% Reconstruct R
measFile = "R2png_coloc_02-Dec-2020_20-56-39";
truthImFile = "truth_ims/R2.png";
fistaOpts.lam1 = 0;
fistaOpts.lam2 = 0;
fistaOpts.maxItr = 30;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;

deconvExp;

%%
clear;

expFname = "R2png_full_03-Dec-2020_00-08-10";
truthImFile = "R2.png";

fistaOpts.lam1 = 2e-2;
fistaOpts.lam2 = 0;
fistaOpts.maxItr = 40;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;

reconExpAbs;

%%

clear;

expFname = "circline2png_full_03-Dec-2020_01-08-12";
truthImFile = "circline2.png";

fistaOpts.lam1 = 1e-2;
fistaOpts.lam2 = 0;
fistaOpts.maxItr = 60;
fistaOpts.tol = 0;
fistaOpts.nonneg = true; 
fistaOpts.showFigs = false;

reconExpAbs;

