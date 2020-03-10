function [dp, ieuler] = emebsddi_dp_multi4b(euler,options)
%EMEBSDDI_DP_MULTI4B Compute dot product by calling EMsoft EMEBSDDI program
% Run EMsoft EMEBSDDI program to simulate EBSD pattern for multiple orientations,
% single set of detector parameters, compute dot product with a multiple
% experimental patterns
% 2/6/20 (Edward Pang, MIT)
%
%%% Inputs: 
% >Np=number of patterns stored in options.img_exp
% >euler=[phi1 PHI phi2
%         phi1 PHI phi2
%         ...          ] (deg)
% >options*
%
% *options is a struct containing the following fields:
% options.homepath = '/home/jonathan/epang/EMsoftfiles/EMdata/';  % path to EMdatapathname
% options.tmppath = 'tmp/';  % path within EMdatapathname for temp files
% options.path = '191212_ZrO2_simdataset/subset1to5/';  % path within EMdatapathname to find .data file
% options.img_exp = 'test.data';    % name of .data file containing experimental patterns (located in path)
% options.ipf_ht = 150;     % height of data set in pattern input file
% options.ipf_wd = 150;     % width
% options.L = 15594.81;    % in microns
% options.xpc = 0.9558;         % in px
% options.ypc = 65.8671;       % in px
% options.thetac = 10.0;   % title angle of the camera (positive below horizontal, degrees)
% options.delta = 59.2;   % CCD pixel size on the scintillator surface (microns)
% options.numsx = 480;    % number of CCD pixels along x and y
% options.numsy = 480;
% options.omega = 0;      % angle between normal of sample and detector
% options.energymin = 10.0;   % energy range in the intensity summation (keV)
% options.energymax = 20.0;
% options.masterfile = 'DItutorial_Ni-master-20kV/Ni-master-20kV.h5';   % master pattern input file, path relative to EMdatapathname
% options.binning = 8;        % binning mode (1, 2, 4, 8)
% options.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
% options.gammavalue = 0.33;  % gamma correction factor
% options.nthreads = 8;
% options.maskpattern = 'y';        % use circular mask? y or n
% options.r = 30;        % radius of circular mask for computing dp
% options.hipassw = 0.05;        % hi pass filter w param; 0.05 is reasonable
% options.nregions = 10;      % # of regions for adaptive histogram equalization
% options.numdictsingle = 1024;   % number of dictionary files arranged in column for dp on GPU (multiples of 16 perform better)
% options.numexptsingle = 1024;   % number of experiment files "
% options.chunk = 2000;     % number of angles to send at once to EMEBSDDI
% options.progress = 1;     % print progress every this many chunks (0 does not show any progress)
%
% Outputs:
% >dp: N_angles x Np array of dp in same order as euler, each col is a map point
% >ieuler: 1 x Np array of indices of euler that correspond to best dp



% Name paths
eulerpath = fullfile(options.homepath,options.tmppath,'euler.txt');
eulerpathshort = fullfile(options.tmppath,'euler.txt');
nmlpath = fullfile(options.homepath,options.tmppath,'EMEBSDDI.nml');
h5path = fullfile(options.homepath,options.tmppath,'DIoutput.h5');
h5pathshort = fullfile(options.tmppath,'DIoutput.h5');
outpathshort = fullfile(options.path,options.img_exp);


% figure out some things
N = options.ipf_wd*options.ipf_ht;  % number of experimental patterns
chunk = options.chunk;
N_angles = size(euler,1);   % how many orientations


% get dp from emebsddi, chunking/padding as necessary
if N_angles <= chunk
    %%% only need to run emebsddi once
    
    if N_angles<chunk
        %%% pad list of angles with dummy orientations and run
        
        % pad list of angles up to chunk
        euler_pad = zeros(chunk,3);  % initialize
        euler_pad(1:N_angles,:) = euler;    % fill in euler angles (padded remain 0)

        % make nml and euler.txt files
        makenml(options,chunk,nmlpath,outpathshort,h5pathshort,eulerpathshort);
        makeeuler(euler_pad,eulerpath);

        % Run EMsoft EMEBSDDI
        [~,~] = system(sprintf('EMEBSDDI %s',nmlpath));

        % read info from .h5 file
        [dp_pad, ieuler] = readh5(h5path);

        % remove dummy values
        dp = dp_pad(1:N_angles,:);
    else
        %%% number of orientations exactly equals chunk size, run list of angles as is

        % make nml and euler.txt files
        makenml(options,N_angles,nmlpath,outpathshort,h5pathshort,eulerpathshort);
        makeeuler(euler,eulerpath);

        % Run EMsoft EMEBSDDI
        [~,~] = system(sprintf('EMEBSDDI %s',nmlpath));

        % read info from .h5 file
        [dp, ieuler] = readh5(h5path);
    end

    % delete input files
    [~,~] = system(sprintf('rm %s',eulerpath));
    [~,~] = system(sprintf('rm %s',nmlpath));
    [~,~] = system(sprintf('rm %s',h5path));
    
else
    %%% need to run emebsddi multiple times and stitch data together
    
    % pad angles if necessary
    if mod(N_angles,chunk)>0
        %%% need to pad list of angles to make a multiple of chunk
        Nchunks = ceil(N_angles/chunk);     % number of chunks to run
        N_angles_pad = Nchunks*chunk;    % total number angles after padding
        euler_pad = zeros(N_angles_pad,3);    % initialize
        euler_pad(1:N_angles,:) = euler;    % fill in euler angles (padded remain 0)
        
        % initialize some arrays
        dp_pad = zeros(N_angles_pad,N);

        % loop through all chunks
        looptime = tic;
        for ii=1:Nchunks
            % make nml and euler.txt files
            makenml(options,chunk,nmlpath,outpathshort,h5pathshort,eulerpathshort);
            makeeuler(euler_pad((ii-1)*chunk+1:ii*chunk,:),eulerpath);

            % Run EMsoft EMEBSDDI
            [~,~] = system(sprintf('EMEBSDDI %s',nmlpath));

            % read info from .h5 file
            [dp_temp, ~] = readh5(h5path);

            % store
            dp_pad((ii-1)*chunk+1:ii*chunk,:) = dp_temp;

            % delete input files
            [~,~] = system(sprintf('rm %s',eulerpath));
            [~,~] = system(sprintf('rm %s',nmlpath));
            [~,~] = system(sprintf('rm %s',h5path));
            
            % print progress
            if mod(ii,options.progress)==0
                percentcomplete = ii/Nchunks;
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
        
        % remove dummy values
        dp = dp_pad(1:N_angles,:);
        
    else
        %%% number of angles is exactly a multiple of chunk, no need to pad
        Nchunks = N_angles/chunk;
        
        % initialize some arrays
        dp = zeros(N_angles,N);

        % loop through all chunks
        looptime = tic;
        for ii=1:Nchunks
            % make nml and euler.txt files
            makenml(options,chunk,nmlpath,outpathshort,h5pathshort,eulerpathshort);
            makeeuler(euler,eulerpath);

            % Run EMsoft EMEBSDDI
            [~,~] = system(sprintf('EMEBSDDI %s',nmlpath));

            % read info from .h5 file
            [dp_temp, ~] = readh5(h5path);

            % store
            dp((ii-1)*chunk+1:ii*chunk,:) = dp_temp;

            % delete input files
            [~,~] = system(sprintf('rm %s',eulerpath));
            [~,~] = system(sprintf('rm %s',nmlpath));
            [~,~] = system(sprintf('rm %s',h5path));
            
            % print progress
            if mod(ii,options.progress)==0
                percentcomplete = ii/Nchunks;
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
    end

    % figure out ieuler
    ieuler = zeros(1,N);    % initialize
    for ii=1:N
        [~,imax] = max(dp(:,ii));  % get index of max dp angle
        ieuler(ii) = imax;
    end
end

end




function makenml(options,N_angles,nmlpath,outpathshort,h5pathshort,eulerpathshort)
% make nml file for emebsddi

fid = fopen(nmlpath,'w');
fprintf(fid,' &EBSDIndexingdata\n');
fprintf(fid,' indexingmode = ''dynamic'',\n');
fprintf(fid,' Notify = ''off'',\n');
fprintf(fid,' ipf_ht = %g,\n',options.ipf_ht);    % height of data set in pattern input file
fprintf(fid,' ipf_wd = %g,\n',options.ipf_wd);          % width
fprintf(fid,' ROI = 0 0 0 0,\n');       % leave all at 0 for full field of view
fprintf(fid,' stepX = 1.0,\n');         % X and Y sampling step sizes
fprintf(fid,' stepY = 1.0,\n');
fprintf(fid,' nnk = %g,\n',N_angles);          % # of top dot products to save (this defines the max # of orientations)
fprintf(fid,' nnav = 1,\n');
fprintf(fid,' nosm = 1,\n');
fprintf(fid,' maskfile = ''undefined'',\n');
fprintf(fid,' maskpattern = ''%s'',\n',options.maskpattern);
fprintf(fid,' maskradius = %.0f,\n',options.r);
fprintf(fid,' hipassw = %g,\n',options.hipassw);
fprintf(fid,' nregions = %g,\n',options.nregions);

fprintf(fid,' ncubochoric = 40,\n');
fprintf(fid,' L = %g,\n',options.L);
fprintf(fid,' thetac = %g,\n',options.thetac);
fprintf(fid,' delta = %g,\n',options.delta);
fprintf(fid,' numsx = %g,\n',options.numsx);
fprintf(fid,' numsy = %g,\n',options.numsy);
fprintf(fid,' xpc = %g,\n',options.xpc);
fprintf(fid,' ypc = %g,\n',options.ypc);
fprintf(fid,' omega = %g,\n',options.omega);
fprintf(fid,' energymin = %g,\n',options.energymin);
fprintf(fid,' energymax = %g,\n',options.energymax);
fprintf(fid,' energyaverage = 0,\n');
fprintf(fid,' spatialaverage = ''n'',\n');
fprintf(fid,' beamcurrent = 150.0,\n');
fprintf(fid,' dwelltime = 100.0,\n');
fprintf(fid,' binning = %g,\n',options.binning);
fprintf(fid,' scalingmode = ''%s'',\n',options.scalingmode);
fprintf(fid,' gammavalue = %g,\n',options.gammavalue);

fprintf(fid,' exptfile = ''%s'',\n',outpathshort);
fprintf(fid,' inputtype = ''Binary'',\n');
fprintf(fid,' HDFstrings = '''' '''' '''' '''' '''' '''' '''' '''' '''' '''',\n');

fprintf(fid,' tmpfile = ''EMEBSDDict_tmp.data'',\n');
fprintf(fid,' keeptmpfile = ''n'',\n');
fprintf(fid,' datafile = ''%s'',\n',h5pathshort);
fprintf(fid,' ctffile = ''undefined'',\n');
fprintf(fid,' avctffile = ''undefined'',\n');
fprintf(fid,' eulerfile = ''%s'',\n',eulerpathshort);

fprintf(fid,' dictfile = ''undefined'',\n');

fprintf(fid,' masterfile = ''%s'',\n',options.masterfile);

fprintf(fid,' numdictsingle = %g,\n',options.numdictsingle);
fprintf(fid,' numexptsingle = %g,\n',options.numexptsingle);
fprintf(fid,' nthreads = %g,\n',options.nthreads);
fprintf(fid,' platid = 1,\n');
fprintf(fid,' devid = 1,\n');
fprintf(fid,' multidevid = 1 0 0 0 0 0 0 0,\n');
fprintf(fid,' usenumd = 1,\n');

fprintf(fid,' /\n');
fclose(fid);

end


function makeeuler(euler,eulerpath)
% Create euler.txt file for emebsddi

fid = fopen(eulerpath,'w');
fprintf(fid,'eu\n');
fprintf(fid,'%.0f\n',size(euler,1));
for ii=1:size(euler,1)
    fprintf(fid,' %g %g %g\n',euler(ii,1),euler(ii,2),euler(ii,3));
end
fclose(fid);

end


function [dp, ieuler] = readh5(h5path)
% read in .h5 file outputted by emebsddi

% sometimes can't find h5path immediately after running EMEBSDDI, these lines of code try to remedy that
t = tic;    
while exist(h5path,'file')==0 && toc(t)<1
    pause(0.1);     % pause 0.1 seconds and try again
end
if exist(h5path,'file')==2
    TopMatchIndices = h5read(h5path,'/Scan 1/EBSD/Data/TopMatchIndices');
    TopDotProductList = h5read(h5path,'/Scan 1/EBSD/Data/TopDotProductList');
    NumExptPatterns = h5read(h5path,'/Scan 1/EBSD/Data/NumExptPatterns');
else
    error('h5path not found!');
end

% Chop empty data
if size(TopMatchIndices,2) > NumExptPatterns
    TopMatchIndices(:,NumExptPatterns+1:end) = [];
    TopDotProductList(:,NumExptPatterns+1:end) = [];
end

% Resort dot products in same order as EulerAngles
dp = zeros(size(TopDotProductList,1),NumExptPatterns);  % initialize
ieuler = TopMatchIndices(1,:);     % only save euler index for best dp
for ii=1:NumExptPatterns
    dp(TopMatchIndices(:,ii),ii) = TopDotProductList(:,ii);
end

end
