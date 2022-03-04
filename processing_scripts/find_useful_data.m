% find_useful_data.m
%
% This routine searches through the meta data (meta data file name 
% to be specified by user) for all instruments that meet the user's 
% specific criteria. Such criteria could be latitude, longitude, 
% depth, or time ranges, type of data included, length or resolution
% of record, etc.
% 
% The selected records are then (optionally) scanned for redundant
% records. The user can specify which of the available databases 
% should be pulled from first, in the event of redundancies.
%
% Finally, all instrument records at a given latitude, longitude, 
% overlapping in time, are grouped into moorings (will be records
% from different depths).
% 
% The selected meta data will be written into a new file.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following section is to be modified by the user and
% specifies file names, data selection criteria, and the priority
% to be applied in the event of redundancies:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gmacmdroot = fullfile(pwd, '../');  
addpath(fullfile(gmacmdroot, 'GMACMD/GMACMDlib'));
addpath(fullfile(gmacmdroot, 'GMACMD/GMACMDlib/tools'));
addpath(fullfile(gmacmdroot, 'GMACMD/GMACMDlib/qc'));

% 1. NAME OF META DATA FILE to search within:
meta_infile = fullfile(gmacmdroot, 'GMACMD/GMACMDdata/metadata_ALL.mat');

% 2. NAME OF OUTPUT META DATA FILE:
meta_outfile = '../data/internal/mooring_metadata.mat';


% 3. SELECTION CRITERIA
% To turn checking of a variable on set check flag to 1, to turn off set
% check flag to 0. E.g. To turn latitude checking on, set 
% criteria.chlat = 1, to turn off set criteria.chlat = 0. When checking of
% a variable is switched on, set assosicated parameters to desired values.

% Search by latitude (degrees north in [-90,90]):
criteria.chlat = 1;
criteria.lat_min = -90;
criteria.lat_max = 90;

% Seach by longitude (degrees east in [0,360)):
criteria.chlon = 1;
criteria.lon_min = 0;
criteria.lon_max = 360;

% Search by time range (see "help datenum" in matlab):
% No records to begin after time_max; no records to end before time_min.
criteria.chtime = 1;
criteria.time_min = datenum(1993,1,1,0,0,0); 
criteria.time_max = datenum(2021,12,31,0,0,0);

% Seach by record length (days):
criteria.chlength = 1;
criteria.reclength_min = 60;

% Search by resolution (minutes) (if criteria.chdelt = 1, then time series
% with irregular time vectors will still be selected if
% criteria.chdelt_irreg = 1, otherwise they will be rejected):
criteria.chdelt = 1;
criteria.chdelt_irreg = 0;
criteria.delt_max = 60;

% Search by instrument depth range (meters, positive below surface,
% negative above): 
criteria.chdepth = 1;
criteria.depth_min = 0;
criteria.depth_max = 200;

% Search by seafloor depth (positive meters) (if chseafloor_nan = 1,
% interpolated sea floor values from bottom topography are used when sea
% floor meta data variable is NaN):
criteria.chseafloor = 1;
criteria.chseafloor_nan = 1;
criteria.seafloor_min = 1000;
criteria.seafloor_max = 30000;

% Search by instrument height above sea floor (meters):
criteria.chheight = 0;
criteria.height_max = 500;
criteria.height_min = 0;

% Desired (required) variables:
% Set a row vector containing natural numbers denoting which
% variables are desired. The numbers are related to variables as,
% 1 = u, 
% 2 = v, 
% 3 = w, 
% 4 = pressure, etc.
% (See 'Var_def' in metadata.)
criteria.wantvars = [1, 2];

% 4. REDUNDANCY PROTOCAL
% Should script check for redundant records?
% (1 for YES, 0 for NO [0 returns all records, regardless of redundancies])
criteria.checkredundancies = 1;

% Set the following flag to 1 in order that, when redundancies are found,
% the record with the greatest duration is kept.
criteria.redund_duration_flag = 1;

% Archive preference. If redund_duration_flag ~= 1, then arhive preference
% is used to select which record of a set of redundant records is kept.
% Further, if redund_duration_flag = 1, and a number of records in a
% redundant set share the largest duration in the set, then archive
% preference is used to select which of those records is kept.
% Let:
% 1 = OSU
% 2 = Wunsch 1997
% 3 = UOPG
% 4 = BODC
% 5 = FAOC
% 6 = Misc
% 7 = NODC
% 8 = IFREMER
% 9 = CSIRO
% 10 = JODC
% 11 = ANIMATE
% 12 = IEO
% 13 = NIOZ
% 14 = ENEA
% Then redundancyopt should be a 1x14 vector containing the numbers
% 1,...,14 in the order that reflects sub-archive preference.
% E.g. redundancyopt = [2,1,3,4,5,6,7,8,9,10,11,12,13,14], means give
% preference to Wunsch 1997, then OSU, then UOPG, etc.
criteria.redundancyopt = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];

% What level of QC do we want the records to have:
% 0 : all records
% 1 : all non-suspicious records (including those not QCed)
% 2 : only QCed records including those flagged as suspicious
% 3 : only QCed records that have not been flagged as suspicious
criteria.qcflag = 3;

% What topograhy file do we wish find_useful_instruments to use? (The
% topography is used to check for instruments on land and to find sea floor
% depths for records that have no recorded sea floor depth (i.e.,
% Seafloor_depth = NaN):
% 1: Smith and Sandwell
% 2: ETOPO 2 min. topography interpolated onto Aviso grid
% 3: Smith and Sandwell interpolated onto Aviso grid
% 4: ETOPO 2 min. topography
% 5: SRTM Plus 30 sec. topograhy (slows down routine significantly)
criteria.topoflag = 2;

% 5. PARAMETERS:
eps.xy = 0.005;      % the tolerance for deciding if two instruments
                     % are at the same location, in degrees E or N (~200m)
eps.t = 0.01;        % time tolerance in days for beginning and end times
eps.z = 2;           % the tolerance for deciding if two instruments
                     % are at the same depth
qc_tol.len = 20;     % section length in days used in quality control
qc_tol.std = 5e-4;   % tolerance for deciding if time series is 
                     % unresonably constant over sections



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the end of the user-modified section. The script below
% should not need to be altered.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: find all the instruments that satisfy the user's particular
% criteria:
find_useful_instruments(meta_infile, meta_outfile, criteria, eps);

% Step 2: further quality control check records
qc_records2(meta_outfile, meta_outfile, criteria.wantvars, qc_tol);

% Step 3: group those instruments into moorings (same lat, lon, time,
% different depths):
find_useful_moorings(meta_outfile, eps);






