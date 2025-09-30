function [Out] = GetYeoNetworks(nP, nNet, varargin)
%GETYEONETWORKS Get Schaefer atlas
%
%   [Out] = GETYEONETWORKS(nP, nNet, varargin)
%
%   Inputs:
%     - nP: number of parcels
%     - nNet: number of networks
%
%   Outputs:
%     - Out: structure with fields:
%       - Explanation: explanation of the fields
%       - Names: (1, nNet) cellstr, network names
%       - Hem: (1, nP) cellstr, hemisphere ('L' or 'R') before sorting
%       - Net: (1, nP) double, network indices before sorting
%       - ParcelNames: (1, nP) cellstr, parcel names before sorting
%       - SortInd: (1, nP) double, indices to sort parcels according to Names
%       - SortNet: (1, nP) double, network indices after sorting with SortInd
%
%   This function simply reads the Schaefer parcellation `info.txt` file
%   under `data/atlas/` folder and extracts the relevant information. Therefore,
%   make sure you have the correct version of the `info.txt` file there and be
%   consistent with the version you used to parcellate your data. Version
%   inconsistency will lead to wrong parcel ordering and naming.

% Input
if nargin < 2
    nNet = 17;
end
assert(nNet == 7 | nNet == 17, 'nNet must be 7 or 17');
if mod(nP, 100)
    nP = 100 * floor(nP / 100);
    warning('Ignoring subcortical regions');
end

% Atlas directory
atlasDir = fullfile('data', 'atlas');

Out = struct();
Out.Explanation = ['ParcelNames, Hem and Net gives the names, hemisphere ' ...
    'and network indices of each parcel before sorting. ' ...
    'SortInd gives the indices to sort parcels according to Names. ' ...
    'SortNet gives the network indices after sorting with SortInd.'];
if nNet == 17
    Out.Names={'VisCent','VisPeri','SomMotA','SomMotB','TempPar',...
        'DorsAttnA','DorsAttnB','SalVentAttnA','SalVentAttnB',...
        'ContA','ContB','ContC','DefaultA','DefaultB','DefaultC',...
        'LimbicA','LimbicB'};
else
    Out.Names={'Vis','SomMot','DorsAttn','SalVentAttn',...
        'Cont','Default','Limbic'};
end

% Full parcel names with index (e.g., 'LH_DorsAttn_Post_1')
fname = fullfile(atlasDir, ...
    sprintf('Schaefer2018_%dParcels_%dNetworks_order_info.txt', nP, nNet));
fileID = fopen(fname);
C = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
C = strrep(C{1}(1:2:end), sprintf('%dNetworks_', nNet), '');
C = C';  % Make row vector
Out.ParcelNames = C;

% Hemisphere and network indices
syllables = cellfun(@(x) strsplit(x, '_'), C, 'UniformOutput', false);
Out.Hem = cellfun(@(x) x{1}(1), syllables, 'UniformOutput', false);  % Strip away "H" in "LH"
Out.Net = cellfun(@(x) find(strcmpi(x{2}, Out.Names)), syllables);

% Sort parcels according to network
[Out.SortNet, Out.SortInd] = sort(Out.Net);  % In a same network, LH comes before RH as in the info file

end