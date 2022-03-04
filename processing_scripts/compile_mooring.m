gmacmdroot = fullfile(pwd, '../');  
% addpath(fullfile(gmacmdroot, 'GMACMD/GMACMDlib'));
% addpath(fullfile(gmacmdroot, 'GMACMD/GMACMDlib/tools'));
% addpath(fullfile(gmacmdroot, 'GMACMD/GMACMDlib/qc'));

meta = load('../data/internal/mooring_metadata.mat');

imooring = 44;
ninst = meta.Mooring_ninst(imooring);

file_number = meta.Mooring_instruments(imooring, 1:ninst);
files = meta.Source_mat(file_number);

I0 = load(fullfile(gmacmdroot, files{1}));

M.time = get_time_vector(I0.begintime, I0.endtime, I0.increment, I0.u);
ntime = length(M.time);

M.depth = NaN(1, ninst);
M.u = NaN(ntime, ninst);
M.v = NaN(ntime, ninst);

for i = 1:numel(files)
    file = files{i};
    IN = load(fullfile(gmacmdroot, file));
    M.u(:, i) = IN.u;
    M.v(:, i) = IN.v;
    M.depth(i) = IN.idepth;
end

if ninst > 1
    figure(101)
    pcolor(M.time, M.depth, M.u')
    shading flat
    datetick
    axis ij
elseif ninst == 1
    figure(101)
    plot(M.time, M.u)
    datetick
end

ncfn = 'test.nc';

nccreate(ncfn, 'lon');
nccreate(ncfn, 'lat');
nccreate(ncfn, 'bdepth');
nccreate(ncfn, 'time', 'Dimensions', {'time' ntime});
nccreate(ncfn, 'depth', 'Dimensions', {'depth' ninst});
nccreate(ncfn, 'u', ...
         'Dimensions', {'time', ntime, 'depth', ninst}, ...
         'FillValue', NaN);
nccreate(ncfn, 'v', ...
         'Dimensions', {'time', ntime, 'depth', ninst}, ...
         'FillValue', NaN);

ncwrite(ncfn, 'bdepth', I0.sfdepth);
ncwriteatt(ncfn, 'bdepth', 'long_name', 'bottom depth')
ncwriteatt(ncfn, 'bdepth', 'standard_name', 'sea_floor_depth_below_sea_surface')
ncwriteatt(ncfn, 'bdepth', 'units', 'm')

ncwrite(ncfn, 'lon', I0.longitude);
ncwriteatt(ncfn, 'lon', 'long_name', 'longitude')
ncwriteatt(ncfn, 'lon', 'standard_name', 'longitude')
ncwriteatt(ncfn, 'lon', 'units', 'degree_east')

ncwrite(ncfn, 'lat', I0.latitude);
ncwriteatt(ncfn, 'lat', 'long_name', 'latitude')
ncwriteatt(ncfn, 'lat', 'standard_name', 'latitude')
ncwriteatt(ncfn, 'lat', 'units', 'degree_north')

ncwrite(ncfn, 'time', posixtime(datetime(M.time, 'ConvertFrom', 'datenum')));
ncwriteatt(ncfn, 'time', 'long_name', 'posix time')
ncwriteatt(ncfn, 'time', 'standard_name', 'time')
ncwriteatt(ncfn, 'time', 'units', 's')

ncwrite(ncfn, 'depth', M.depth);
ncwriteatt(ncfn, 'depth', 'long_name', 'depth')
ncwriteatt(ncfn, 'depth', 'standard_name', 'depth')
ncwriteatt(ncfn, 'depth', 'units', 'm')

ncwrite(ncfn, 'u', M.u);
ncwriteatt(ncfn, 'u', 'long_name', 'eastward velocity')
ncwriteatt(ncfn, 'u', 'standard_name', 'eastward_sea_water_velocity')
ncwriteatt(ncfn, 'u', 'units', 'm s-1')

ncwrite(ncfn, 'v', M.v);
ncwriteatt(ncfn, 'v', 'long_name', 'northward velocity')
ncwriteatt(ncfn, 'v', 'standard_name', 'northward_sea_water_velocity')
ncwriteatt(ncfn, 'v', 'units', 'm s-1')

% % Check records
% endtimes = zeros(meta.Mooring_ninst(imooring), 1);
% for i = 1:numel(files)
%     file = files{i};
%     load(fullfile(gmacmdroot, file), 'endtime');
%     endtimes(i) = endtime;
% end
% 
% fprintf('STD of end time = %3.1f\n', std(endtimes))

% function M = compile_mooring(gmacmdroot, files)
% 
% 
% end

% Test pairing

% gmacmdroot = fullfile(pwd, '../'); 
% meta = load('../data/internal/mooring_metadata.mat');
% 
% vars = zeros(size(meta.Var_def));
% vars([1, 2]) = 1;
% tol.lat = 0.005;
% tol.lon = 0.005;
% tol.idepth = 0.5;
% tol.begintime = 1;   % in days
% tol.endtime = 1;
% pairmeta('EW_global/data/internal/mooring_metadata.mat', 'mooring_metadata_paired.mat', vars, tol)