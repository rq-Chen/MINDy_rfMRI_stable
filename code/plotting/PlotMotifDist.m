function f = PlotMotifDist(motifs, centroids)
%PLOTMOTIFDIST Plot the distribution of each parcel's activation in motifs
if nargin < 2
    centroids = mean(motifs, 3)';
end

nParcels = size(motifs, 1);
nFigPerMotif = nParcels / 100;
nMotifs = size(motifs, 2);
atlas = GetYeoNetworks(nParcels);
SortInd = atlas.SortInd;
parcelNames = atlas.ParcelNames;
nNetworks = numel(atlas.Names);
c = distinguishable_colors(nNetworks);
f = cell(nMotifs, nFigPerMotif);
for i = 1:nMotifs
    for iFig = 1:nFigPerMotif
        ff = figure();
        f{i, iFig} = ff;
        t = tiledlayout(8, 13);
        tid = [1:89 92:102];  % Tile ID
        for j = 1:100
            ax = nexttile(tid(j));
            hold on;
            jj = SortInd((iFig - 1) * 100 + j);
            cc = c(atlas.Net(jj), :);
            histogram(motifs(jj, i, :), 'Normalization', 'pdf', 'FaceColor', cc);
            xline(mean(motifs(jj, i, :)), '--k', 'LineWidth', 2);
            xlim([-1 1] * max(abs(xlim)));
            ax.XColor = cc;
            ax.YColor = cc;
            title(parcelNames{jj}, 'Color', cc, 'Interpreter', 'none', 'FontSize', 7);
        end
        nexttile(90, [2 2]);
        PlotYeoSurface(centroids(i, :));
        colormap turbo
        colorbar;
        title('Centroid');

        xlabel(t, 'Parcel activity');
        ylabel(t, 'Probability density function')
        if nFigPerMotif > 1
            title(t, sprintf('Marginal distribution of motif %d (%d of %d)', ...
                i, iFig, nFigPerMotif));
        else
            title(t, sprintf('Marginal distribution of motif %d', i));
        end
        ff.WindowState = 'maximized';
    end
end

end

