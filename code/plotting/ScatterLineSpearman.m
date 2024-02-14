function[Out,r]=ScatterLineSpearman(X,Y,varargin)
ttTemp=and(and(~isnan(X(:)),~isnan(Y(:))),and(~isinf(X(:)),~isinf(Y(:))));
Out=scatter(X(ttTemp),Y(ttTemp),varargin{:});
lsline
r=nancorr(X(:),Y(:),'type','Spearman');
title(strcat('\rho=',num2str(r)));
grid on
end