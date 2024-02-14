%% Infinite Period Bifurcation model

clear; clc; close all;

scSize = 100;
scColor = 'g';


%% Simulate the limit cycle

T = 200;

allMu = 1.1:-0.01:1.01;
Nmu = length(allMu);
allCoeffs = nan(2, 2, Nmu);
for i = 1:Nmu
    sims = ipbMdl(allMu(i), [1; rand() * 2 * pi], T);
    sims = [cos(sims(2, :)); sin(sims(2, :))];
    allCoeffs(:, :, i) = pca(sims');
end
allCoeffs = allCoeffs .* sign(allCoeffs(2, 1, :));

figure('Units', 'inches', 'Position', [0 0.5 18 8]);
gscatter(squeeze(allCoeffs(1, 1, :)), squeeze(allCoeffs(2, 1, :)), 1:Nmu, parula(Nmu));


%% Vector field plot

allx = -1.2:0.1:1.2;
ally = allx';
[allx, ally] = meshgrid(allx, ally);
allr = sqrt(allx .^ 2 + ally .^ 2);
allthetasin = ally ./ allr;
allthetacos = allx ./ allr;

drawArrow = @(x, y, varargin) quiver(x(1), y(1), x(2)-x(1), y(2)-y(1), 0, varargin{:});

allMu = [1.1 1.01 0.99 0.9];
Nmu = length(allMu);

t = tiledlayout(2, Nmu, "TileSpacing", "tight");
for i = 1:Nmu
    mu = allMu(i);
    allrdot = allr .* (1 - allr .^ 2);
    allthetadot = -abs(allthetasin) + mu;
    allxdot = allrdot .* allthetacos - allr .* allthetasin .* allthetadot;
    allydot = allrdot .* allthetasin + allr .* allthetacos .* allthetadot;

    nexttile(i)
    hold on;
    streamslice(allx, ally, allxdot, allydot);
    if mu > 0 && mu < 1
        scatter([sqrt(1 - mu ^ 2) -sqrt(1 - mu ^ 2)], [mu -mu], scSize, 'r', 'filled');
        scatter([-sqrt(1 - mu ^ 2) sqrt(1 - mu ^ 2)], [mu -mu], scSize, 'y', 'filled');
        scatter(nan, nan, scSize, scColor, "filled");
    elseif mu > 1
        sims = ipbMdl(allMu(i), [1; rand() * 2 * pi], 100, 0.05);
        sims = [cos(sims(2, :)); sin(sims(2, :))];
        scatter(sims(1, :), sims(2, :), [], 'k', 'filled');        
        scatter([0; 0], [1; -1], scSize, scColor, "filled");
    end
    axis equal
    axis tight
    title(['\mu = ' sprintf('%.2f', mu)]);
    subtitle('Trajectories')

    nexttile(i + Nmu);
    hold on;
    quiver(allx, ally, allxdot, allydot, 5);
    if mu > 0 && mu < 1
        s1 = scatter([sqrt(1 - mu ^ 2) -sqrt(1 - mu ^ 2)], [mu -mu], scSize, 'r', 'filled');
        s2 = scatter([-sqrt(1 - mu ^ 2) sqrt(1 - mu ^ 2)], [mu -mu], scSize, 'y', 'filled');
        s3 = scatter(nan, nan, scSize, scColor, "filled");
        s4 = scatter(nan, nan, [], 'k', 'filled');
    elseif mu > 1
        sims = ipbMdl(allMu(i), [1; rand() * 2 * pi], 100, 0.05);
        sims = [cos(sims(2, :)); sin(sims(2, :))];
        scatter(sims(1, :), sims(2, :), [], 'k', 'filled');        
        scatter([0; 0], [1; -1], scSize, scColor, "filled");
    end
    axis equal
    axis tight
    subtitle('Vector field')
end
lgd = legend([s1 s2 s3 s4], {'Stable fixed point', 'Saddle node', 'Ghost', 'Simulated samples'});
lgd.Layout.Tile = 'east';
fontsize(gcf(), 16, 'points');
title(t, "r' = r(1 - r^2), \theta' = \mu - abs(sin \theta)", 'FontWeight', 'bold', ...
    'FontSize', 20);
PrintAsSeen(fullfile('figures', 'SNIC'), '-dpng', '-r300');


%% Model function
function sims = ipbMdl(mu, x, T, dt)
%IPBMDL Simulation an infite period bifurcation model
%
%   The model is in polar axis: $r' = r(1 - r^2), \theta' = \mu -
%   abs(\theta)$
%
%   Inputs:
%   - mu: a scaler around 1
%   - x0: (2, N) matrix. x0(1, :) should be $r$ and x0(2, :) is $\theta$
%   - T: simulation time
%   - dt: deltaT, default is 0.1
%
%   Output:
%   - sims: (2, T / dt + 1, N) matrix

    if nargin < 4
        dt = 0.1;
    end
    N = size(x, 2);
    L = T / dt + 1;
    sims = nan(2, L, N);
    for l = 1:L
        sims(:, l, :) = reshape(x, 2, 1, N);
        x(1, :) = x(1, :) + (x(1, :) .* (1 - x(1, :) .^ 2)) * dt;
        x(2, :) = x(2, :) + mu - abs(sin(x(2, :)));
    end
end