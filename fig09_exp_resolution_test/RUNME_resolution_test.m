clear; close all; clc;

addpath("../lib");
addpath("../export_fig");
savePath = sprintf("../../figs/fig09_exp_resolution_test");

testSep = 0.5;
resolution_test;

testSep = 1.0;
resolution_test;

testSep = 2.0;
resolution_test;