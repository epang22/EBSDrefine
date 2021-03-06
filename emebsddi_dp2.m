function dp = emebsddi_dp2(L,xpc,ypc,euler,options)
% EMEBSDDI_DP2
% Run EMsoft EMEBSDDI program to compute dot products for multiple orientations,
% single set of detector parameters, compute dot product with a single experimental pattern
% 1/29/20 (Edward Pang, MIT)
%
%%% Inputs: 
% >L (um)
% >xpc (px)
% >ypc (px)
% >euler=[phi1 PHI phi2
%         phi1 PHI phi2
%         ...          ] (deg)
% >options*
%
% *options is a struct containing the following fields:
% options.homepath = '/home/jonathan/epang/EMsoftfiles/EMdata/';  % path to EMdatapathname
% options.tmppath = 'tmp/';  % path within EMdatapathname for temp files and to find img_exp
% options.img_exp = 'test.data';    % name of .data file containing experimental pattern
% options.thetac = 5.0;   % title angle of the camera (positive below horizontal, degrees)
% options.delta = 50.0;   % CCD pixel size on the scintillator surface (microns)
% options.numsx = 480;    % number of CCD pixels along x and y
% options.numsy = 480;
% options.omega = 0;      % angle between normal of sample and detector
% options.energymin = 10.0;   % energy range in the intensity summation (keV)
% options.energymax = 25.0;
% options.masterfile = '180713_13p5Ce_25kV_master/ZrO2-13p5CeO2_tetragonal-master-25kV.h5';   % master pattern input file, path relative to EMdatapathname
% options.binning = 1;        % binning mode (1, 2, 4, 8)
% options.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
% options.gammavalue = 0.33;  % gamma correction factor
% options.maskpattern = 'y';        % use circular mask? y or n
% options.r = 230;        % radius of circular mask for computing dp
% options.nthreads = 16;
% options.numdictsingle = 1024;   % number of dictionary files arranged in column for dp on GPU (multiples of 16 perform better)
% options.numexptsingle = 1024;   % number of experiment files "
% options.nnk = 1000;          % # of top dot products to save (this defines the max # of orientations)



% Name paths
eulerpath = strcat(options.homepath,options.tmppath,'euler.txt');
eulerpathshort = strcat(options.tmppath,'euler.txt');
nmlpath = strcat(options.homepath,options.tmppath,'EMEBSDDI.nml');
h5path = strcat(options.homepath,options.tmppath,'DIoutput.h5');
h5pathshort = strcat(options.tmppath,'DIoutput.h5');
outpathshort = strcat(options.tmppath,options.img_exp);


N_angles = size(euler,1);   % Figure out how many angles




% Create EMEBSDDI.nml file
fid = fopen(nmlpath,'w');
fprintf(fid,' &EBSDIndexingdata\n');
fprintf(fid,' indexingmode = ''dynamic'',\n');
fprintf(fid,' Notify = ''off'',\n');
fprintf(fid,' ipf_ht = 1,\n');    % height of data set in pattern input file
fprintf(fid,' ipf_wd = 1,\n');          % width
fprintf(fid,' ROI = 0 0 0 0,\n');       % leave all at 0 for full field of view
fprintf(fid,' stepX = 1.0,\n');         % X and Y sampling step sizes
fprintf(fid,' stepY = 1.0,\n');
fprintf(fid,' nnk = %g,\n',options.nnk);
fprintf(fid,' nnav = 1,\n');
fprintf(fid,' nosm = 1,\n');
fprintf(fid,' maskfile = ''undefined'',\n');
fprintf(fid,' maskpattern = ''%s'',\n',options.maskpattern);
fprintf(fid,' maskradius = %.0f,\n',options.r);
fprintf(fid,' hipassw = 0.05,\n');      % hi pass filter w param; 0.05 is reasonable
fprintf(fid,' nregions = 10,\n');       % # of regions for adaptive histogram equalization

fprintf(fid,' ncubochoric = 40,\n');
fprintf(fid,' L = %g,\n',L);
fprintf(fid,' thetac = %g,\n',options.thetac);
fprintf(fid,' delta = %g,\n',options.delta);
fprintf(fid,' numsx = %g,\n',options.numsx);
fprintf(fid,' numsy = %g,\n',options.numsy);
fprintf(fid,' xpc = %g,\n',xpc);
fprintf(fid,' ypc = %g,\n',ypc);
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


    
% Create euler.txt file
fid = fopen(eulerpath,'w');
fprintf(fid,'eu\n');
fprintf(fid,'%.0f\n',size(euler,1));
for ii=1:N_angles
    fprintf(fid,' %g %g %g\n',euler(ii,1),euler(ii,2),euler(ii,3));
end
fclose(fid);

% Run EMsoft EMEBSDDI
[~,~] = system(sprintf('EMEBSDDI %s',nmlpath));
% system(sprintf('EMEBSDDI %s',nmlpath));     % verbose


% Read in data
% sometimes can't find h5path immediately after running EMEBSDDI, these lines of code try to remedy that
t = tic;    
while exist(h5path,'file')==0 && toc(t)<5
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
% if size(TopMatchIndices,1) > N_angles
%     TopMatchIndices(N_angles+1:end) = [];
%     TopDotProductList(N_angles+1:end) = [];
% end

% Resort dot products in same order as EulerAngles
dp = zeros(length(TopDotProductList),1);  % initialize
for ll = 1:length(TopDotProductList)
    index = TopMatchIndices(ll);
    dp(index) = TopDotProductList(ll);
end
    


% delete input files
[~,~] = system(sprintf('rm %s',eulerpath));
[~,~] = system(sprintf('rm %s',nmlpath));
[~,~] = system(sprintf('rm %s',h5path));

