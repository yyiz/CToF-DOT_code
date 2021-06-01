clear; close all; clc;

addpath("../lib");
addpath("../export_fig");

datdir = "../dat/12_9_20_vis_shift_psf_Jacobian";
savedir = "../../figs/fig04_confocal_vis_jacobian";
savename = "shift_invar";
savenameY_plot = "../../figs/fig04_confocal_vis_jacobian/shift_invar_Ysec";
savenameX_plot = "../../figs/fig04_confocal_vis_jacobian/shift_invar_Xsec";
fname = "J";

finds = 0:3;

f1 = figure();
f2 = figure('Position', [114.5000 256 539.5000 276.5000]);
f3 = figure('Position', [114.5000 256 539.5000 276.5000]);
for i = 1:length(finds)
    if (i == 1)
        plotRowX = 16;
        plotRowY = 16;
        plotColor = 'r';
    elseif (i == 2)
        plotRowX = 16;
        plotRowY = 9;
        plotColor = 'g';
    elseif (i==3)
        plotRowX = 11;
        plotRowY = 27;
        plotColor = 'b';
    else
        plotRowX = 23;
        plotRowY = 11;
        plotColor = 'k';
    end
    
    load(sprintf("%s/%s%d.mat", datdir, fname, finds(i)));
    visJ = reshape(sum(J, 1), Jheaders.VOX_W, Jheaders.VOX_L, Jheaders.VOX_H);
    visJ = flipud(visJ'); % Flipd in Y for visualization, X and Y dimensions are flipped
    figure(f1)
    subplot(2,2,i);
    imagesc(visJ);
    yline(plotRowY, plotColor, 'LineWidth', 3);
    xline(plotRowX, plotColor, 'LineWidth', 3);
    
    srcdetX = (Jheaders.SRC_ORIG(1) - Jheaders.VOX_ORIGX)/Jheaders.VOX_SIDELEN;
    srcdetY = (-Jheaders.SRC_ORIG(2) - Jheaders.VOX_ORIGY)/Jheaders.VOX_SIDELEN; % Down is up in programming world
    hold on;
    
    set(gca, "XTickLabel", []); set(gca, "YTickLabel", []);
    pbaspect([1 1 1]);
    
    
    middleInd = Jheaders.VOX_L/2 - 1;
    xAx = 1:32;
    figure(f2);
    hold on;
    psfIY = visJ(plotRowY,:);
    [~, maxInds] = maxk(psfIY,2);
    maxInd = mean(maxInds);
    shft = middleInd - maxInd;
    
    newAx = xAx - shft;
    psfIY = interp1(xAx, psfIY, newAx, 'nearest', 'extrap');
    
    plot(xAx, psfIY, plotColor, 'LineWidth', 3);
    xlim([1, 32]);
    set(gca, 'FontSize', 16);
    set(gca, 'FontName', 'Times New Roman');
    
    figure(f3);
    hold on;
    psfIX = visJ(:,plotRowX);
    [~, maxInds] = maxk(psfIX,2);
    maxInd = mean(maxInds);
    shft = middleInd - maxInd;
    
    newAx = xAx - shft;
    psfIX = interp1(xAx, psfIX, newAx);

    plot(xAx, psfIX, plotColor, 'LineWidth', 3);
    xlim([1, 32]);
    set(gca, 'FontSize', 16);
    set(gca, 'FontName', 'Times New Roman');
end

figname = sprintf("%s/%s.fig", savedir, savename);
pngname = sprintf("%s/%s.png", savedir, savename);
fignameY = sprintf("%s/%s.fig", savedir, savenameY_plot);
pngnameY = sprintf("%s/%s.png", savedir, savenameY_plot);
fignameX = sprintf("%s/%s.fig", savedir, savenameX_plot);
pngnameX = sprintf("%s/%s.png", savedir, savenameX_plot);

if ~exist(savedir, 'dir')
    mkdir(savedir);
end

if (exist(figname, 'file'))
    fprintf("File exists, overwrite?\n");
    keyboard;
end
savefig(f1, figname);
export_fig(f1, pngname, '-m3', '-transparent', '-png');
savefig(f2, fignameY);
export_fig(f2, pngnameY, '-m3', '-transparent', '-png');
savefig(f3, fignameX);
export_fig(f3, pngnameX, '-m3', '-transparent', '-png');

