clear; close all; clc;

defaultPath = "../dat";
savedir = "../../figs/suppl_optical_param_measure";
addpath("../lib");
addpath("../export_fig");

saveName = "mus_tpsf_mat";
mFnum = 5:12;
mFpath = sprintf("%s/3_31_21_tpsf_broadening/tempfiles", defaultPath);
outPath = sprintf("%s/%s", savedir, saveName);

if ~exist(outPath, 'dir')
    mkdir(outPath);
end

mPath = sprintf("%s/background_tpsf%d.txt", mFpath, mFnum(1));
[mHeaders, nMHeaderLines] = parseHeader(mPath);

mus_vec = zeros(length(mFnum),1);
tpsf_mat = zeros(length(mFnum), mHeaders.NBINS);

DTIME = (mHeaders.TIME_MAX - mHeaders.TIME_MIN)/mHeaders.NBINS;
timeAx = (mHeaders.TIME_MIN+(DTIME/2)):DTIME:mHeaders.TIME_MAX;

legendVec = [];

f = figure(); hold on;
for i = 1:length(mFnum)
    mPath = sprintf("%s/background_tpsf%d.txt", mFpath, mFnum(i));
    [mHeaders, nMHeaderLines] = parseHeader(mPath);
    mus_vec(i) = mHeaders.SCAT_VEC;
    tpsf_mat(i,:) = readmatrix(mPath, "NumHeaderLines", nMHeaderLines);
    plot(tpsf_mat(i,:), 'LineWidth', 2);
    legendVec = [legendVec sprintf("\\mu_s=%d mm^{-1}", mus_vec(i))];
end

legend(legendVec);
set(gca, 'FontSize', 16);
set(gca, 'FontName', 'Times New Roman');
set(gca,'YTickLabel',[]);

xlabel('Time (ps)')

savefig(f, outPath);
% saveas(f, sprintf("%s.png", outPath));
export_fig(f, outPath, '-nocrop', '-transparent', '-m3', '-png');

