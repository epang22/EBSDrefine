%RUNEBSDREFINE_CHECKANGLES
% Use this script to check list of compiled euler angles before running 
% EBSDrefine to figure out what values for 'eulerminappear','angletol',
% 'fitmax', and 'reduceeuler' to use
% Fill in INPUT PARAMETERS section with desired parameters
% 2/20/20 (Edward Pang, MIT)


clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.path = 'testdata_EBSDrefine/';  % path within EMdatapathname to find .ang file
data.inputfile = 'Scan3_Mod.ang';     % .ang file name of Hough data (located in path)

data.phaseid = 2;    % which phase to consider (as numbered in ang file)

% refined grid settings for interpolation
data.angletol = 1.2;   % find all unique euler angles at least this many degrees apart as grid centers
data.eulerminappear = 5;     % minimum number of times euler angle must appear in Hough data to be considered
data.fitmax = 0.5;  % maximum fit value for that euler angle to be considered
data.reduceeuler = 1;     % eliminate symmetry related euler angles? 1=yes, 0=no
    % =1 is slower to find unique euler angles, but saves time later for indexing
    % Set to 1 unless the finding euler angles step takes a really long time (>30min)

% Define pseudosymmetry variants to check (comment out if you don't want to check pseudosym)
data.pseudosym = [
    1 1 0 90;
    1 -1 0 90
    ];     % [axis_x axis_y axis_z angle(deg)]



%%% Parameters you don't need to change often %%%
% Paths for this computer
data.homepath = '/home/jonathan/epang/EMsoftfiles/EMdata/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DO NOT EDIT BELOW THIS LINE %%%


EBSDrefine_checkangles(data)

