function [Out, Dat, dDat, X1, dX1] = MINDy_from_Deconv(Dat, mdlType, derivType, nBatch, svCP, Pre, ParStr)
%MINDY_FROM_DECONV A wrapper for MINDy_Base or MINDy_Linear for comparison
%
%   [Out, Dat, dDat] = MINDy_from_Deconv(Dat, mdlType)
%
%   Inputs:
%     Dat: (nParcels, nTR) matrix or (1, nSess) cell array of such matrices,
%       the DECONVOLED, PREPROCESSED and ZSCORED time series with NO NaNs.
%     mdlType: 'Base' or 'Linear' or 'NoDiag'
%     derivType: 'two-point' (default) or 'one-point' or 'center'
%     nBatch: number of batches for training (default 20000)
%     svCP: whether to save checkpoint (default true)
%     Pre & ParStr: hyperparameter structures from ChosenPARSTR
%
%   Outputs:
%     Out: MINDy output structure with optional parameter update history
%     Dat: (1, nSess) cell, discarding the last two frames in each cell
%     dDat: the derivative corresponding to Dat
%     X1: (1, 1) cell containing the good frames from Dat across all sessions
%     dX1: the derivative corresponding to X1
%
%   In order to understanda the influence of early stopping, here the models
%   will be trained for much longer time. The parameter update history will
%   be saved in |RecW|, |RecA| and |RecD| fields of the output. The mean
%   squared error (for each parcel in each batch) will be saved in |E|.

% Inputs
if nargin < 6
    ChosenPARSTR;
end
if nargin < 5
    svCP = true;
end
if nargin < 4
    nBatch = 20000;  % Use much longer training time (originally 5000)
end
if nargin < 3
    derivType = 'two-point';
end
if nargin < 2
    mdlType = 'Base';
end
if ~iscell(Dat)
    Dat={Dat};
end
nX=size(Dat{1},1);

% Parameters
ParStr.BatchSz=300;
ParStr.NBatch=nBatch;
if svCP  % Record parameter update history (every 50 batches)
    ParStr.RecW='y';
    ParStr.RecA='y';
    ParStr.RecDdiff='y';
else
    ParStr.RecW='n';
    ParStr.RecA='n';
    ParStr.RecDdiff='n';
end
ParStr.L2SpPlsEN=0;  % No L2 normalization on the low rank component
% Rescale dimensions of low-rank component proportionately to the number of parcels
ParStr.wPC=ceil(ParStr.wPC*nX/419);
% HRF and TR stuff
ParStr.H1min=5;ParStr.H1max=7;ParStr.H2min=.7;ParStr.H2max=1.3;ParStr.H1Rate=.1;ParStr.H2Rate=.1;
Pre.TR=.72;
% Slope for linear activation function
slope = 0.62;  % Match the nonlinearity across randn(): 0.62; match around 0: 1.3

% Derivatives
switch derivType
    case 'two-point'
        dDat=cellfun(@(xx)(convn(xx,[1 0 -1]/2,'valid')),Dat,'UniformOutput',0);
        Dat=cellfun(@(xx)(xx(:,1:end-2)),Dat,'UniformOutput',0);
    case 'one-point'
        dDat=cellfun(@(xx)(convn(xx,[1 -1],'valid')),Dat,'UniformOutput',0);
        Dat=cellfun(@(xx)(xx(:,1:end-1)),Dat,'UniformOutput',0);
    case 'center'
        dDat=cellfun(@(xx)(convn(xx,[1 0 -1]/2,'valid')),Dat,'UniformOutput',0);
        Dat=cellfun(@(xx)(xx(:,2:end-1)),Dat,'UniformOutput',0);
end


if strcmpi(mdlType, 'Linear')
    Out=MINDy_Linear(Dat,dDat,Pre,ParStr,[],slope);
elseif strcmpi(mdlType, 'Base')
    Out=MINDy_Base(Dat,dDat,Pre,ParStr);
elseif strcmpi(mdlType, 'NoDiag')
    Out=MINDy_Base(Dat,dDat,Pre,ParStr,logical(ones(nX)-eye(nX)));
else
    error('Unknown model type: %s', mdlType);
end

% Global regression on W and D to get rid of global regularization bias
if strcmpi(mdlType, 'Linear')  % Add a function handle to compute derivative for MINDy_Inflate
    Out=MakeMINDyFunction_Linear(Out);
else
    Out=MakeMINDyFunction(Out);
end
[X1,dX1]=MINDy_Censor(Out,Dat,dDat);  % Use good frames for inflation and correlation
Out=MINDy_Inflate(Out,X1,dX1,'n');  % No need for robust regression

% Recalculate goodness-of-fit
Out.Corr=DiagCorr(Out.FastFun([X1{:}])',[dX1{:}]');
% Out=rmfield(Out,{'FastFun','Tran','dTran'});  % Function handle generates warnings
Out.Pre=Pre;Out.ParStr=ParStr;

end