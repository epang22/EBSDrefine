function EBSDrefine_nonpseudo( data )
%EBSDREFINE_NONPSEUDO
% Read in Hough data from .ang file, create refined grid
% around each angle and pseudosym variants, feed into EMEBSDDI, extract
% data and interpolate refined orientations, output refined .ang file
% Only works for non-pseudosymmetric materials
% 2/20/20 (Edward Pang, MIT)
%
%%% Input: 'data', a struct containing the following fields:
% data.homepath = '/home/jonathan/epang/EMsoftfiles/EMdata/';  % path to EMdatapathname
% data.tmppath = 'tmp/';   % path within EMdatapathname for creation of temp files
% 
% data.path = '180627ELP33scan3/fullmap/';  % path within EMdatapathname to find .ang file and to save output
% data.inputfile = 'Scan3_Mod.ang';     % .ang file name of Hough data (located in path)
% data.img_exp = 'Scan3.data';         % .data file containing all patterns (located within path)
% data.data.binning = 2;        % binning mode (1, 2, 4, 8). Simulated image size will be numsx/binning x numsy/binning
% 
% data.h5output = '200207_ebsd_zro2pseudosym_onephase3e_monoclinic_angletol1p2_eulerminappear5_fitmax1_da2_Na4.h5';   % Name of .h5 output file (located in path)
% data.angoutput = '200207_ebsd_zro2pseudosym_onephase3e_monoclinic_angletol1p2_eulerminappear5_fitmax1_da2_Na4.ang';   % Name of .ang output file (located in path)
% 
% data.phaseid = 1;    % which phase to consider (as numbered in ang file)
% 
% % refined grid settings for interpolation
% data.angletol = 1.2;   % find all unique euler angles at least this many degrees apart as grid centers
% data.eulerminappear = 5;    % minimum number of times euler angle must appear in Hough data to be considered
% data.fitmax = 1;  % maximum fit value for that euler angle to be considered
% data.reduceeuler = 1;     % eliminate symmetry related euler angles? 1=yes, 0=no
% 
% data.delta_angles = 2;   % check +/- this many degrees for refined orientation grids
% data.N_angles = 4;   % this many steps from orientation to orientation+delta_angles (2N+1 points in each dim)
% 
% 
% %%% Info about data set
% % pattern center
% data.L = 15733.10;    % in microns (snobfit monoclinic)
% data.xpc = 0.9415;         % in px
% data.ypc = 92.8340;       % in px
% 
% 
% %%% EMsoft parameters for simulating patterns and computing dot products
% data.thetac = 5;   % tilt angle of the camera (positive below horizontal, degrees)
% data.delta = 50;   % CCD pixel size on the scintillator surface (microns)
% data.numsx = 480;    % number of CCD pixels along x and y
% data.numsy = 480;
% data.omega = 0;      % angle between normal of sample and detector
% data.energymin = 10.0;   % energy range in the intensity summation (keV)
% data.energymax = 25.0;
% data.masterfile = '180713_13p5Ce_25kV_master/ZrO2-13p5CeO2_monoclinic-master-25kV.h5';   % master pattern input file, path relative to EMdatapathname
% data.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
% data.gammavalue = 0.33;  % gamma correction factor
% data.maskpattern = 'y';        % use circular mask? y or n
% data.hipassw = 0.05;    % hi pass filter w param; 0.05 is reasonable
% data.nregions = 10;     % # of regions for adaptive histogram equalization
% data.r = min(data.numsx,data.numsy)/(2*data.binning);        % radius of circular mask (after binning)
% data.nthreads = 7; % for EMEBSDDI
% 
% % number of dictionary/experiment files arranged in column for dp on GPU (multiples of 16 perform better)
% data.numdictsingle = 192;
% data.numexptsingle = 192;
% data.chunk = 2000;      % max number of orientations to send to emebsddi at once
% data.progress = 1;      % print progress every this many chunks
% 
% % How many iterations to print progress
% data.progressinterp = 2000;   % attempts to print progress each map point of interpolation
% data.verbose = 0;    % print warnings to screen? 1=yes, 0=no



%%% Extract input parameters from struct 'data'
homepath = data.homepath;
path = data.path;
tmppath = data.tmppath;
inputfile = data.inputfile;
h5output = data.h5output;
angoutput = data.angoutput;
phaseid = data.phaseid;
angletol = data.angletol;
delta_angles = data.delta_angles;
N_angles = data.N_angles;
fitmax = data.fitmax;
eulerminappear = data.eulerminappear;
reduceeuler = data.reduceeuler;
progressinterp = data.progressinterp;
verbose = data.verbose;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


totaltime = tic;

% clear tmp folder
[~,~] = system(sprintf('rm %s%s*',homepath,tmppath));

% Error checking
h5path = strcat(homepath,path,h5output);    % full path to h5 output file
if exist(h5path,'file')==2
    error('File already exists! Rename .h5 output file.');
end
angpath = strcat(homepath,path,angoutput);    % full path to ang output file
if exist(angpath,'file')==2
    error('File already exists! Rename .ang output file.');
end


%%% Read in data
disp('Reading in data...');

% read in ang file data
inputpath = fullfile(homepath,path,inputfile);   % full path of .ang file to read in
[euler,x,y,IQ,~,phase,fit,phaseinfo,grid,PChough] = loadang(inputpath,1);

% get OIM symmetry ID of phase of interest
for ii=1:size(phaseinfo,1)
    if phaseinfo{ii,1}==phaseid
        phaseindex = ii;
        symid = phaseinfo{ii,3};   % get OIM symmetry ID
        break   % found it, don't need to keep checking
    end
end

% store some info
ipf_ht = grid{6};
ipf_wd = grid{4};
data.ipf_ht = ipf_ht;
data.ipf_wd = ipf_wd;
if strcmp(grid{1},'HexGrid')
    issquare = 0;
else
    issquare = 1;
end
data.issquare = issquare;

% figure out indices of dummy patterns
if issquare==0
    idummy = 2*ipf_wd:2*ipf_wd:ipf_wd*ipf_ht;
end



%%% Loop through Euler angles and compile list of unique Euler angles
disp('Finding unique Euler angles...');

eulerlist1 = [];    % keep track of all unique euler angles
eulertally = [];    % keep track of how many times euler angle appears (within tol)

N = size(euler,1);     % number of map points

% loop through each map point
for ii=1:N
    if phase(ii)==phaseid   % only for phase of interest
        if fit(ii)<fitmax   % consider this map point only if fit < this threshold
            euler1 = euler(ii,1:3);   % extract euler angles for this map point

            % compare to existing angles in list
            counter = 1;    % if counter=1 at end of eulerlist, then add euler angle to list
            for jj=1:size(eulerlist1,1)
                euler2 = eulerlist1(jj,:);   % extract euler angles for this iteration
                misorient = misorientation(euler1,euler2);

                if misorient < angletol
                    counter = 0;    % don't add eulerangle to list
                    eulertally(jj) = eulertally(jj) + 1;
                    break   % already found a match, don't need to check rest of eulerlist
                end
            end

            % if not within angletol of any euler angle in list, add euler angles of this map point to list
            if counter==1
                eulerlist1(end+1,:) = euler1;
                eulertally(end+1) = 1;
            end
        end
    end

    % print progress
    if mod(ii,1000)==0
        fprintf('  %g of %g completed...\n',ii,N);
    end
end



% remove symmetry related euler angles, if specified
if reduceeuler==1
    disp('Removing symmetry related Euler angles...');

    N2 = size(eulerlist1,1);    % number of angles in eulerlist1 (picked out by misorientation)
    eulerlist1_disorient = [];
    eulertally_disorient = [];

    % loop through eulerlist1
    for ii=1:N2
        euler1 = eulerlist1(ii,:);   % extract euler angle for this iter
        eulertally_temp = eulertally(ii);   % extract number of map points for this euler angle

        % compare to existing angles in eulerlist1_disorient
        counter = 1;    % if counter=1 at end, then add euler angle to list     
        for jj=1:size(eulerlist1_disorient,1)
            euler2 = eulerlist1_disorient(jj,:);   % extract euler angles for this iteration
            [~,disorient] = disorientation(euler1,euler2,symid);

            if disorient < angletol
                counter = 0;    % don't add eulerangle to list
                eulertally_disorient(jj) = eulertally_disorient(jj) + eulertally_temp;
                break   % already found a match, don't need to check rest of eulerlist1_disorient
            end
        end

        % if not within angletol of any euler angle in eulerlist1_disorient, add to list
        if counter==1
            eulerlist1_disorient(end+1,:) = euler1;
            eulertally_disorient(end+1) = eulertally_temp;
        end

        % print progress
        if mod(ii,100)==0
            fprintf('  %g of %g completed...\n',ii,N2);
        end
    end
else
    % copy to new variables
    eulerlist1_disorient = eulerlist1;
    eulertally_disorient = eulertally;
end



% remove angles in eulerlist that appear less than eulerminappear times
eulerlist = eulerlist1_disorient(eulertally_disorient>=eulerminappear,:);
if size(eulerlist,1)<1
    error('Zero unique euler angles found. Reduce eulerminappear.');
end
Neuler = size(eulerlist,1);     % number of grid centers




%%% add grid of angles to refine
gridsize = (2*N_angles+1)^3;    % number of points in each grid
eulerlist_grids = zeros(Neuler*gridsize,3);   % contains all refined grids about each angle in eulerlist
for ii=1:Neuler
    euler1 = eulerlist(ii,:);
    [test_angles, test_angles_1Drodrigues] = createorigrid(delta_angles, N_angles, euler1);   % Create test grid of angles
    eulerlist_grids(gridsize*(ii-1)+1:ii*gridsize,:) = test_angles;
end
eulerlist_grids_all = eulerlist_grids;  % compile dictionary: original variant only




%%% Dictionary index to find best euler angle for each map point
disp('Performing refined dictionary indexing...');

[dp_all, ieuler] = emebsddi_dp_multi4b(eulerlist_grids_all, data);

% if data hexagonal, remove data from dummy patterns
if issquare==0
    ieuler(idummy) = [];
    dp_all(:,idummy) = [];
end



%%% extract data, figure out which grids, interpolate
disp('Performing interpolation...');

% initialize some arrays
%  best interpolated/actual euler angles+dp for each map point and pseudosym variant
%   dim1: euler angles (deg) then dp, dim2: map points
maxinterppoint_all = zeros(4,N);

% warnings. store patternid and pseudosym var# for bad interpolations
warnings = [];  % col1: patternid, col2: var# (original=0), col3: exitflag

% loop through each map point
looptime = tic;     % start timer for loop
for ii=1:N
    % get indexed orientation
    euler_indexed = eulerlist_grids_all(ieuler(ii),:);
    
    % figure out which grid center is closest to indexed orientation
    misorient = 180;    % initialize, keep track of lowest misorientation btw indexed orientation and grid centers
    for jj=1:Neuler
        % check original vars
        misorient_temp = misorientation(euler_indexed,eulerlist(jj,:));
        if misorient_temp<misorient
            igridcenter = jj;
            misorient = misorient_temp;
        end
    end
        
    % orientation of best grid center
    euler1 = eulerlist(igridcenter,:);    % grid center of original variant
    q1 = eu2qu(euler1*pi/180);     % Quaternion of grid center
    
    % interpolate
    dp_grid = dp_all(gridsize*(igridcenter-1)+1:igridcenter*gridsize,ii);     % extract dp for original variant grid
    dp_3d = reshapeto3d(dp_grid, 2*N_angles+1);  % reshape array to 3d
    [maxinterppoint,maxgridpoint,~,~,exitflag] = dpinterp(dp_3d, test_angles_1Drodrigues, q1);

    % check if there are any issues in interpolated data
    if exitflag>0
        if verbose==1
            fprintf('%g. Interpolation warning. Exitflag=%g.\n',ii,exitflag);
        end
        warnings(end+1,:) = [ii 0 exitflag];

        % store best grid point instead of interpolated point
        maxinterppoint_all(:,ii,1) = maxgridpoint;
    else
        % store interpolated point
        maxinterppoint_all(:,ii,1) = maxinterppoint;
    end
    
    
    %%% print progress
    if mod(ii,progressinterp)==0    
        percentcomplete = ii/length(ieuler);
        esttimeremain = toc(looptime)*(1-percentcomplete)/percentcomplete;     % estimated time remaining in seconds

        if esttimeremain<60
            fprintf('    %.1f%% completed. Estimated time remaining: %.1f seconds...\n',100*percentcomplete,esttimeremain);
        elseif esttimeremain<3600
            fprintf('    %.1f%% completed. Estimated time remaining: %.1f minutes...\n',100*percentcomplete,esttimeremain/60);
        else
            fprintf('    %.1f%% completed. Estimated time remaining: %.1f hours...\n',100*percentcomplete,esttimeremain/3600);
        end 
    end
end



% calculate some things
eulerbest = maxinterppoint_all(1:3,:)';     % euler angles of best variant (rad)
dpbest = maxinterppoint_all(4,:)';        % dp of best variant
CI = dpbest - mean(dp_all(:));      % use peak height above bkg as CI
    % almost all comparisons are rubbish comparisons, so should give
    % reasonable estimate of the background noise dp value




% mis-/dis-orientation btw Hough and DI+interp for each map point
misorient_hough = zeros(N,1);   % initialize
disorient_hough = zeros(N,1);   % initialize
for ii=1:N
    [miso, diso] = disorientation(euler(ii,:), eulerbest(ii,:), symid);
    misorient_hough(ii) = miso;
    disorient_hough(ii) = diso;
end

% plot histograms
figure; histogram(disorient_hough); xlabel('Disorientation (deg)'); ylabel('Count');



% calculate some things
fit = 1-dpbest;     % measure of OIM's fit for DI
euler = euler*pi/180;   % convert to rad for storing
eulerbest = eulerbest*pi/180;   % convert to rad for storing
xstar = data.xpc/data.numsx + 0.5;
ystar = data.ypc/data.numsy + 0.5;
zstar = data.L/(data.numsx*data.delta);
phase2 = zeros(N,1);    % phase=0, like OIM ang file for 1 phase



% summarize warnings to screen
fprintf('\n%g warnings.\n',size(warnings,1));



executiontime = toc(totaltime);



% store data
data.euler = euler;
data.eulerlist_all = eulerlist1;    % before filtering by disorientation and eulerminappear
data.eulerlist = eulerlist;
data.eulerlist_grids_all = eulerlist_grids_all;     
data.ieuler = ieuler;
data.dp_all = dp_all;
data.test_angles_1Drodrigues = test_angles_1Drodrigues;
data.maxinterppoint_all = maxinterppoint_all;
data.eulerbest = eulerbest';
data.dpbest = dpbest';
data.CI = CI';
data.fit = fit';
data.info = 'Indexed by ebsd_zro2pseudosym_onephase3e (employing EMsoft EMEBSDDI)';
if ~isempty(PChough)
    data.PC = [xstar ystar zstar PChough(4)];   % DI pattern center, WD from microscope ang
else
    data.PC = [xstar ystar zstar 10];   % use WD=10 as dummy if no data available
end
data.IQ = IQ';     % from Hough
data.phase = phase2';
data.x = x';
data.y = y';

phaseinfo_export{1} = phaseinfo{phaseindex,1};
phaseinfo_export{2} = phaseinfo{phaseindex,2};
phaseinfo_export{3} = phaseinfo{phaseindex,3};
phaseinfo_export{4} = phaseinfo{phaseindex,4};
data.phaseinfo = phaseinfo_export;

if ~isempty(grid)
    data.grid = grid;   % only load to feed into makeang.m if there is actually grid info
end
data.warnings = warnings';
data.executiontime = executiontime;
data.misorient_hough = misorient_hough';
data.disorient_hough = disorient_hough';


% make ang file
makeang(angpath, data, 1);
fprintf('Data stored in ang file: %s\n',angpath);

% save data to hdf5 file
makehdf5_3d(h5path, data);
fprintf('Data stored in h5 file: %s\n\n',h5path);


fprintf('Elapsed time is %.6f seconds.\n',executiontime);




