gmacmdroot = fullfile(pwd, '../');  
addpath(fullfile(gmacmdroot, 'GMACMD/GMACMDlib'));
addpath(fullfile(gmacmdroot, 'GMACMD/GMACMDlib/tools'));
addpath(fullfile(gmacmdroot, 'GMACMD/GMACMDlib/qc'));

meta = load('../data/internal/mooring_metadata.mat');

meanKE = NaN*meta.Mooring_instruments;

% Now loop through the moorings:
for i = 1:size(meta.Mooring_instruments,1)
   
    % Loop trough instruments on a mooring:
    for j = 1:meta.Mooring_ninst(i)
        
        % The instrument number:
        n = meta.Mooring_instruments(i,j);
        
        % Load in the time series for this current meter:
        M = load([gmacmdroot, meta.Source_mat{n}]);
        
        % Now that the above line has loaded u and v into the workspace we
        % can take the mean kinetic energy, noting that there may be NaNs
        % in the time series:
        meanKE(i,j) = 0.5*mean(M.u.^2 + M.v.^2, 'omitnan');
    
    end
    
end

% Now let's make a histogram of the mean kinectic energy below the 90th
% percentile:
figure();
q = ~isnan(meanKE);
x = meanKE(q);
y = prctile(x, 90);
figure(1), clf, set(gca, 'FontSize', 12)
histogram(x(x<y));
xlabel('Mean Kinetic Energy [m^2 s^{-2}]')
ylabel('Frequency')
title('Mean KE below the 90th Percentile')

topo = load(fullfile(gmacmdroot, 'GMACMD/GMACMDtopo/ETOPO2v2c_f4_AvisoGrid.mat'));
figure(2);
contour(topo.NbLongitudes, topo.NbLatitudes, topo.altitude_aviso, [0, 0])
hold on
plot(meta.Mooring_lon, meta.Mooring_lat, 'o')