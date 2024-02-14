function Out = MINDyInt_00(ooP, Start, DnSamp, dt, dW, tEnd, u)
%
% Author: Matthew Singh (updated by Ruiqi Chen)
%
% Simulate data from MINDy models.
%
% Inputs:
%
%   ooP: model structure
%   Start: initial conditions, (nX, 1) for single sim or (nX, k) for k parallel sims
%   DnSamp: downsampling factor for saving data (recommended = 1)
%   dt: time-step size in units TR (recommended = 1)
%   dW: noise std (iid). Shape should be compatible with (nX, tEnd / dt, k).
%     - dW is the noise added at each step
%     - It does NOT autoscale based upon dt -- you'll have to do that manually
%   tEnd: sim length, should be dividable by dt
%   u: optional input vector with size compatible to (nX, tEnd / dt, k). Note: if the
%     model input is "Au" you need to pre-multiply u by A before passing it in.
%
% Output:
%
%   Out: (n, tEnd / dt / DnSamp + 1) for a single sim or (n, tEnd / dt / DnSamp + 1, k) for multiple
%

if nargin < 7 || isempty(u) || ~isnumeric(u)
    u = 0;
end

if ~isempty(ooP.Param{5})
    W = ooP.Param{5};
else  % Only sparse but no low-rank component
    W = ooP.Param{1};
end

A = ooP.Param{2};
D = ooP.Param{6};
B = ooP.Param{3};
A2 = A .^ 2;

% Multiply W by b while dividing the activation phi(x) by b
A2 = A2 ./ (B(:, 1) .^ 2);
B5P = (B(:, 2) + .5) ./ B(:, 1);
B5N = (B(:, 2) - .5) ./ B(:, 1);
% Factor out the dt terms and add xt
DDT = 1 - (dt * D);
WDT = dt * W .* (B(:, 1)');

if ~isempty(ooP.Param{4})
    Cdt = ooP.Param{4} * dt;
else
    Cdt = 0;
end

nX = size(Start, 1);
nDat = size(Start, 2);

tVec = 0:(dt * DnSamp):tEnd;
nT = numel(tVec);

Out = nan(nX, nT, nDat);
Out(:, 1, :) = Start;


%% Simulation

if any(dW, "all")
    noise = dW .* randn(nX, tEnd / dt, nDat);
else
    noise = zeros(nX, tEnd / dt, nDat);
end
u = u + noise;

for i = 1:(tEnd / dt)
    Start = WDT * (sqrt(A2 + (Start + B5P) .^ 2) - sqrt(A2 + (Start + B5N) .^ 2)) ...
        + DDT .* Start + squeeze(u(:, i, :)) + Cdt;
    if mod(i, DnSamp) == 0
        Out(:, 1 + (i / DnSamp), :) = Start;
    end
end

end
