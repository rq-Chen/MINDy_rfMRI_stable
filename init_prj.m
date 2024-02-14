function oldpath = init_prj(other_folders)
%INIT_PRJ Set up paths to scripts and dependecies
%
%   This script will set up the project root `prj_root`.
%   It will add the script folders & subfolders to MATLAB search path,
%   and change the working directory to `prj_root`.
%
%   `other_folders` is an optional input. You can use
%   it to add MINDy modeling package folder into the path.
%   (Note that you should always have only one MINDy package
%   in your path in case of conflicts.)

% Folders to add to path
folders = {'code'; 'cifti-matlab-master'};
if nargin > 0 && ~isempty(other_folders)
    folders = [folders; reshape(cellstr(other_folders), [], 1)];
end

% prj_root = 'path/to/project/root'
[prj_root, ~, ~] = fileparts(mfilename("fullpath"));  % the parent dir of this script
fprintf("Project root: %s\n", prj_root);

oldpath = path();
for i = 1:length(folders)
    fullfolder = fullfile(prj_root, folders{i});
    fprintf("Adding %s and all subdirs to path...\n", fullfolder);
    addpath(genpath(fullfolder))
end

% Change directory
fprintf("Changing the working directory to project root...\n")
cd(prj_root)

end