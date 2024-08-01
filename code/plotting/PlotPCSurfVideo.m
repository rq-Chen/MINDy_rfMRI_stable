function [] = PlotPCSurfVideo(mdl, fname, dW, t, showFC)
%PLOTPCSURFGIF Plot the evolution of states in PC1/2 space and on the Yeo surface
%
%   Inputs:
%     mdl: MINDy structure
%     fname: filename to save video to
%     dW: standard deviation of noise to add to states
%     t: time point to plot, 1D index array
%     showFC: whether to show time-varying functional connectivity


%% Parameters

% For replicability
rng(0);
useNoiseFreePC = true;

% Defaults & Constants
if nargin < 5 || isempty(showFC)
    showFC = false;
end
if nargin < 4 || isempty(t)
    t = 1:5:300;
end
if nargin < 3 || isempty(dW)
    dW = sqrt(1 - mdl.Corr .^ 2);
end
if nargin < 2 || isempty(fname)
    fname = 'PCSurf.mp4';
end
nParcels = size(mdl.Param{5}, 1);
nSims = 30; % Number of simulations to run for the distribution
T = max(t);  % Simulation length

fr = 6;  % Frame rate
plotSimID = 3;  % Which simulation to plot
nt = length(t);  % Number of time points to plot

tvFCwin = 20;  % Window size for time-varying functional connectivity (in TRs)


%% Simulation

% Get attractors and PCs
sims = MINDyInt_00(mdl, randn(nParcels, nSims), 1, 1, dW, T - 1);  % (nParcels, T, nSims)
sims = reshape(sims, nParcels, []);
[eq, lc, nfSims] = GetEqLC(mdl);
if ~useNoiseFreePC
    % mu = mean(sims, 2);  % (nParcels, 1)
    mu = zeros(nParcels, 1);
    [coeff, score] = pca(sims', 'NumComponents', 2, 'Centered', false);  % (nParcels, 2), (T * nSims, 2)
else
    nfSims = reshape(nfSims, nParcels, []);
    % mu = mean(nfSims, 2);  % (nParcels, 1)
    mu = zeros(nParcels, 1);
    coeff = pca(nfSims', 'NumComponents', 2, 'Centered', false);
    score = (sims - mu)' * coeff;
end

% Keep the simulation to plot
fullscore = score;  % Full distribution
sims = sims(:, 1 + T * (plotSimID - 1):T * plotSimID);  % (nParcels, T)
score = score(1 + T * (plotSimID - 1):T * plotSimID, :);  % (T, 2)

% Project attractors to PC space
if ~isempty(eq)
    eqScore = (eq - mu)' * coeff;  % (nEq, 2)
end
if size(lc, 2)
    lcScore = cellfun(@(x) (x - mu)' * coeff, lc, 'UniformOutput', false);  % (1, nLC) cell of (nPeriod, 2)
end

% Get time-varying functional connectivity (windowed correlation matrix)
if showFC
    tvFC = nan(nParcels, nParcels, nt);
    for i = 1:nt
        tvFC(:, :, i) = corr(sims(:, t(i):t(i) + tvFCwin - 1)');
    end
end

% Keep the frames to plot
sims = sims(:, t);  % (nParcels, nt)
score = score(t, :);  % (nt, 2)


%% Draw the first frame

f = figure('Units', 'inches', 'Position', [1 1 14 6], 'Color', 'w');
tlo = tiledlayout(1, 2 + showFC, 'Padding', 'compact');

% Distribution of states
nexttile;
hold on;
scatter(fullscore(:, 1), fullscore(:, 2), 10, 'k', 'filled', 'o', ...
    'MarkerFaceAlpha', 0.05, 'MarkerEdgeAlpha', 0.05);

% Axes
% xline(0, 'k--');
% yline(0, 'k--');

% Limit sets
if ~isempty(eq)
    scatter(eqScore(:, 1), eqScore(:, 2), 250, [238 102 119] / 256, 'filled', 'pentagram');
end
if ~isempty(lc)
    for i = 1:numel(lc)
        scatter(lcScore{i}(:, 1), lcScore{i}(:, 2), 10, [204 187 68] / 256, 'filled', 'o');
    end
end

% Current state
sc = scatter(score(1, 1), score(1, 2), 200, [34 136 51] / 256, 'filled', 'hexagram');

% Labeling
xlabel('PC1');
ylabel('PC2');
axis tight
axis equal;
xticks([]);
yticks([]);
axis off;
title('Phase space', 'FontSize', 16)
subtitle(' ')

% Surface plot
nexttile;
cl = quantile(sims(:), [0.01 0.99]);
PlotYeoSurface(sims(:, 1));
clim(cl);
drawnow;
% Required due to some weird bug in the lighting
cla;
[~, mmpL, mmpR, sL, sR] = PlotYeoSurface(sims(:, 1));
clim(cl);
% end of weird bug fix
title('Cortical activation', 'FontSize', 16)
subtitle(' ')

% Time-varying functional connectivity
if showFC
    nexttile;
    MyNetMatYeo(tvFC(:, :, 1));
    nm = gca().Children(end);
    cl2 = quantile(tvFC(:), [0.01 0.99]);
    % cl2 = clim();
    clim(cl2);
    colormap parula;
    % cb = colorbar;
    % cb.Layout.Tile = 'east';
    axis equal
    xlim([0.5 nParcels + 0.5]);
    ylim([0.5 nParcels + 0.5]);
    title('Time-varying functional connectivity', 'FontSize', 16)
    subtitle(' ')
    % Wierd bug fix
    nm.CData = tvFC(:, :, 1);
    drawnow
end

% Title
title(tlo, sprintf('TR %d', t(1)), 'FontWeight', 'bold', 'FontSize', 18);

% Save the first frame
drawnow;
im = frame2im(getframe(f));
if exist(fname, 'file')
    warning([fname ' exists and will be overwritten!'])
    delete(fname);
end
v = VideoWriter(fname, 'MPEG-4');
v.FrameRate = fr;
open(v);
writeVideo(v, im);


%% Draw the rest of the frames

for i = 2:nt
    % Update the current state
    sc.XData = score(i, 1);
    sc.YData = score(i, 2);

    % Update the surface plot
    nexttile(2);
    cla;
    PlotYeoSurface(sims(:, i), mmpL, mmpR, sL, sR);
    clim(cl);

    % Update the time-varying functional connectivity
    if showFC
        nexttile(3);
        nm.CData = tvFC(:, :, i);
    end

    % Update title
    title(tlo, sprintf('TR %d', t(i)), 'FontWeight', 'bold', 'FontSize', 18);
    
    % Save the frame
    drawnow;
    im = frame2im(getframe(f));
    writeVideo(v, im);
end
close(f);
close(v);


end

