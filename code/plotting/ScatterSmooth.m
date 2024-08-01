function [sc] = ScatterSmooth(Xv, Yv, varargin)
%SCATTERSMOOTH  Scatter plot with regression line and confidence interval

sc = scatter(Xv, Yv, varargin{:});
hold on;

mdl = fitlm(Xv, Yv);
xl = xlim();
x = linspace(xl(1), xl(2), 100)';
[y, yci] = predict(mdl, x);
plot(x, y, 'k');

ylow = yci(:, 1);
yhigh = yci(:, 2);
p = fill([x; flip(x)], [ylow; flip(yhigh)], 'k', 'FaceAlpha', 0.3);
p.EdgeColor = 'none';
% text(x(5), ylow(5), sprintf("R^2 = %.2f", mdl.Rsquared.Ordinary), "FontSize", 14, ...
%     "HorizontalAlignment", "left", "VerticalAlignment", "top", 'Margin', 0.1, ...
%     "Interpreter", "tex", "BackgroundColor", [1 1 1 0.5]);
text(x(95), yhigh(95), sprintf("R^2 = %.2f", mdl.Rsquared.Ordinary), "FontSize", 14, ...
    "HorizontalAlignment", "right", "VerticalAlignment", "bottom", 'Margin', 0.1, ...
    "Interpreter", "tex", "BackgroundColor", [1 1 1 0.5]);
fprintf("p = %e\n", mdl.Coefficients{2, "pValue"});

end