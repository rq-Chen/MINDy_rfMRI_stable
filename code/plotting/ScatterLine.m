function[Out,r]=ScatterLine(X,Y,varargin)
%% Scatterplot and correlation for two variables of arbitrary size (vectorized)
%% varargin takes additional plotting options for scatter


ttTemp=and(and(~isnan(X(:)),~isnan(Y(:))),and(~isinf(X(:)),~isinf(Y(:))));
Out=scatter(X(ttTemp),Y(ttTemp),varargin{:});
lsline
if isinteger(X)
    X=single(X);
end
if isinteger(Y)
    Y=single(Y);
end
r=(nancorr(X(:),Y(:)));
title(strcat('r=',num2str(r)));
grid on
end