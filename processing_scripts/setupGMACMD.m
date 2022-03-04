% setupGMACMD.m
%
% Sets the variable, gmacmdroot, and adds GMACMDlib/ and GMACMDlib/tools/
% to path.

gmacmdroot = '../';  

path(path, fullfile(gmacmdroot, 'GMACMD/GMACMDlib'));
path(path, fullfile(gmacmdroot, 'GMACMD/GMACMDlib/tools'));
path(path, fullfile(gmacmdroot, 'GMACMD/GMACMDlib/qc'));
