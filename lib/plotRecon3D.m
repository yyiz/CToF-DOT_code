function f = plotRecon3D(scene3D, varargin)

name_list = {'VOX_ORIG', 'VOXDIM', 'edgeAlpha', 'cbar',...
    'xlab', 'ylab', 'zlab', 'rgb', 'normalizeIm'};
optargs = {[0, 0, 0], [1, 1, 1], 0, true, 'x (mm)', 'y (mm)','z (mm)', [0 255, 0], true};
optargs = parseVarargin(varargin, name_list, optargs);
[VOX_ORIG, VOXDIM, edgeAlpha, cbar, xlab, ylab, zlab, rgb, normalizeIm] = optargs{:};
[VOX_L, VOX_W, VOX_H] = size(scene3D);

if (normalizeIm)
    maxVal = max(scene3D(:));
    scene3D = scene3D ./ maxVal;
end

f = figure();
h1 = axes;
set(h1, 'Zdir', 'reverse');
axis equal;
view([0, 0]);
xlabel(xlab); ylabel(ylab); zlabel(zlab);

rax = linspace(1, rgb(1)./255, 84)';
gax = linspace(1, rgb(2)./255, 84)';
bax = linspace(1, rgb(3)./255, 84)';
colormap([rax, gax, bax]);


% Display heatmap
xBounds = [VOX_ORIG(1), VOX_ORIG(1) + VOX_W*VOXDIM(1)];
yBounds = [VOX_ORIG(2), VOX_ORIG(2) + VOX_L*VOXDIM(2)];
zBounds = [VOX_ORIG(3), VOX_ORIG(3) + VOX_H*VOXDIM(3)];

vol3d('CData', scene3D, 'XData', xBounds, 'YData', yBounds, 'ZData', zBounds, 'EdgeAlpha', edgeAlpha);
numTicks = 11;
if (cbar)
    colorbar('Ticks', linspace(0,1,numTicks), 'TickLabels', linspace(0,maxVal,numTicks));
    caxis([0 1]);
end

end