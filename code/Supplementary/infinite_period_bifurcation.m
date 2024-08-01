%% Infinite Period Bifurcation model

clear; clc; close all;

scSize = 100;
scColor = 'g';
nSims = 30;
T = 100;
dt = 0.05;


%% Vector field plot

allx = -1.2:0.1:1.2;
ally = allx';
[allx, ally] = meshgrid(allx, ally);
allr = sqrt(allx .^ 2 + ally .^ 2);
allthetasin = ally ./ allr;
allthetacos = allx ./ allr;

drawArrow = @(x, y, varargin) quiver(x(1), y(1), x(2)-x(1), y(2)-y(1), 0, varargin{:});

allMu = [5 1.01 0.99 0.2];
Nmu = length(allMu);

f = figure('Units', 'inches', 'Position', [0 0.5 20 10.0208]);
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
        sims = ipbMdl(allMu(i), [ones(1, nSims); (1:nSims) * 2 * pi / nSims + 0.01], T, dt);
        sims = reshape(sims, 2, []);
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
        sims = ipbMdl(allMu(i), [ones(1, nSims); (1:nSims) * 2 * pi / nSims + 0.01], T, dt);
        sims = reshape(sims, 2, []);
        sims = [cos(sims(2, :)); sin(sims(2, :))];
        scatter(sims(1, :), sims(2, :), [], 'k', 'filled');        
        scatter([0; 0], [1; -1], scSize, scColor, "filled");
    end
    axis equal
    axis tight
    subtitle('Vector field')
end
lgd = legend([s1 s2 s3 s4], {'Stable fixed point', 'Saddle node', 'Slowest point', 'Simulated samples'});
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
%   - sims: (2, T + 1, N) matrix

    if nargin < 4
        dt = 0.1;
    end
    step_per_T = round(1 / dt);
    N = size(x, 2);
    L = T * step_per_T + 1;
    sims = nan(2, T + 1, N);
    for l = 1:L
        if mod(l, step_per_T) == 1
            sims(:, (l - 1) / step_per_T + 1, :) = reshape(x, 2, 1, N);
        end
        x(1, :) = x(1, :) + (x(1, :) .* (1 - x(1, :) .^ 2)) * dt;
        x(2, :) = x(2, :) + mu - abs(sin(x(2, :)));
    end
end