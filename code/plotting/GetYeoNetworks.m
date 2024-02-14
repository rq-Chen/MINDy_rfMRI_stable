function[Out]=GetYeoNetworks(nP,varargin)
%% After pre-calculating with d082018preCalcYeo

%% SortNet gives the network labels AFTER sorting
%% SortInd gives the indices for sorting

%% nP=number of parcels
%% Ordered Networks=order used in Yeo figures

%% Out is a structure containing names, network labels 
%% (with and w/out laterality), sorting indices

aa=fullfile('data', 'atlas');

if ~isempty(varargin)&&strcmpi(varargin{1},'y')
    doShort='y';
else
    doShort='n';
end

gg=load(fullfile(aa,strcat('YeoNet_',num2str(nP),'_',doShort,'.mat')));
Out=gg.gg;
Out.Explanation='SortNet gives indices AFTER sorting';

% Full parcel names with index
fname = fullfile('data', 'atlas', ...
    sprintf('Schaefer2018_%dParcels_17Networks_order_info.txt', nP));
fileID = fopen(fname);
C = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
C = strrep(C{1}(1:2:end), '17Networks_', '');
Out.ParcelNames = C;

end