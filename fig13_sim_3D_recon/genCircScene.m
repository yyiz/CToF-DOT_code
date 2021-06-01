clear; close all; clc;

addpath('../lib');

savedir = "3_21_21_sim_3D_Jacobian";
loaddir = "circ-line-ims";
savename = "truth_circ_32x32x6";
savepath = sprintf("../dat/%s/%s", savedir, savename);

thresh = 0.5;

imX = 500; imY = 500; imZ = 6;
imX_final = 32; imY_final = 32;
im_vol = zeros(imX, imY, imZ);

for i = 1:imZ
    im_i = double(rgb2gray(imread(sprintf("../dat/%s/%s/layer%d.png", savedir, loaddir, i))));
    im_vol(:,:,i) = 1 - (im_i ./ max(im_i(:)));
end

truth = imresize3(im_vol, [imX_final, imY_final, imZ]);
truthInds = truth < thresh;
truth(truthInds) = 0;
truth(~truthInds) = 1;
ftarg = plotRecon3D(truth, 'rgb', [255, 0, 0]);

save(savepath, 'truth');