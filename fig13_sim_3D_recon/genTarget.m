clear; close all; clc;

load("../dat/truth_ims/branch.mat");
savepath = "../dat/3_21_21_sim_3D_Jacobian/truth_32x32x6.mat";

temp = double(OUTPUTgrid);
temp = temp ./ max(temp(:));

temp = imresize3(temp, [32, 32, 6]);

threshInds = temp > 0.01;
temp(threshInds) = 1.0;
temp(~threshInds) = 0.0;

truth = zeros(size(temp));
for i = size(temp,3):-1:1
    truth(:,:,i) = temp(:,:,size(temp,3)+1-i);
end

f_truthslices = figure('Position', [1264 579.5000 560 420]);
for i = 1:size(truth, 3)
    subplot(4,3,i)
    imagesc(truth(:,:,i));
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
end

save(savepath, 'truth');
