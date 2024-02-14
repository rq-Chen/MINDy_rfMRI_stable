function J = GetMINDyJacobian(mdl, X)
%GETMINDYJACOBIAN Calculate the Jacobian of the model |mdl| at |X|
%   Inputs:
%       |mdl|: MINDy structure
%       |X|: (nParcel, 1) or (1, nParcel) vector
    if isrow(X); X = X'; end
    assert(isvector(X), 'Input must be a vector, not matrix!');
    W = mdl.Param{5};
    D = mdl.Param{6};
    AA = mdl.Param{2} .^ 2;
    b = mdl.Param{3}(1);
    bxp = b .* X + .5;
    bxm = bxp - 1;
    ppsipx = b .* (bxp ./ sqrt(AA + bxp .^ 2) - (bxm ./ sqrt(AA + bxm .^ 2)));
    J = W .* (ppsipx') - diag(D);
end