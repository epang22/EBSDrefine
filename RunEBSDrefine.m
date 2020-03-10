%RUNEBSDREFINE
% Use this script to perform EBSD data refinement
% Fill in INPUT PARAMETERS section with desired parameters
% 2/20/20 (Edward Pang, MIT)

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define path for input and output data
data.path = 'testdata_EBSDrefine/';  % path within EMdatapathname to find .ang file and to save output

% Input data
data.inputfile = 'Scan3_Mod.ang';     % .ang file name of Hough data (located in path)
data.img_exp = 'Scan3.data';         % .data file containing all patterns (located within path)
data.binning = 2;        % binning mode (1, 2, 4, 8). Simulated image size will be numsx/binning x numsy/binning
data.phaseid = 2;    % which phase to consider (as numbered in ang file header)

% Output data
data.h5output = '200223_EBSDrefine_test_tetragonal.h5';   % Name of .h5 output file (located in path)
data.angoutput = '200223_EBSDrefine_test_tetragonal.ang';   % Name of .ang output file (located in path)

% refined grid settings for interpolation
data.angletol = 1.2;   % find all unique euler angles at least this many degrees apart as grid centers
data.eulerminappear = 5;    % minimum number of times euler angle must appear in Hough data to be considered
data.fitmax = 0.5;  % maximum fit value for that euler angle to be considered
data.reduceeuler = 1;     % eliminate symmetry related euler angles? 1=yes, 0=no

data.delta_angles = 2;   % check +/- this many degrees for refined orientation grids
data.N_angles = 2;   % this many steps from orientation to orientation+delta_angles (2N+1 points in each dim)

% pattern center
data.L = 15733.10;    % in microns
data.xpc = 0.9415;         % in px
data.ypc = 92.8340;       % in px

% Define pseudosymmetry variants to check (comment out if you don't want to check any pseudosymmetry variants)
data.pseudosym = [
    1 1 0 90;
    1 -1 0 90
    ];     % [axis_x axis_y axis_z angle(deg)]

% EMsoft parameters for the master pattern
data.energymin = 10.0;   % energy range in the intensity summation (keV)
data.energymax = 25.0;
data.masterfile = '180713_13p5Ce_25kV_master/ZrO2-13p5CeO2_tetragonal-master-25kV.h5';   % master pattern input file, path relative to EMdatapathname

% EMsoft EMEBSDDI parameters
% *You need to check that these parameters work using TestErrorsTimings.m
% *For whatever reason, EMsoft gives an error for some combinations.
data.nthreads = 7; % for EMEBSDDI
data.numdictsingle = 192;   % # of dictionary files arranged in column for dp on GPU (multiples of 16 perform better)
data.numexptsingle = 192;   % # of experiment files arranged in column for dp on GPU (multiples of 16 perform better)
data.chunk = 2000;      % max number of orientations to send to EMEBSDDI at once

% How many iterations to print progress
data.progress = 1;      % print progress every this many chunks
data.progressinterp = 2000;   % print progress every this many map points during interpolation
data.verbose = 0;    % print warnings to screen? 1=yes, 0=no (they will be saved no matter what)



%%% Parameters you don't need to change often %%%
% Paths for this computer
data.homepath = '/home/jonathan/epang/EMsoftfiles/EMdata/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)
data.tmppath = 'tmp/';   % path within EMdatapathname for creation of temp files

% Detector parameters
data.thetac = 5;   % tilt angle of the camera (positive below horizontal, degrees)
data.delta = 50;   % CCD pixel size on the scintillator surface (microns)
data.numsx = 480;    % number of CCD pixels along x and y
data.numsy = 480;
data.omega = 0;      % angle between normal of sample and detector

% Some EMsoft parameters for dot product computations
data.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
data.gammavalue = 0.33;  % gamma correction factor
data.maskpattern = 'y';        % use circular mask? y or n
data.hipassw = 0.05;    % hi pass filter w param; 0.05 is reasonable
data.nregions = 10;     % # of regions for adaptive histogram equalization
data.r = min(data.numsx,data.numsy)/(2*data.binning);        % radius of circular mask (after binning)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %



if isfield(data,'pseudosym')
    EBSDrefine_pseudo(data)
else
    EBSDrefine_nonpseudo(data)
end




