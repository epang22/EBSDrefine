%RUNEBSDREFINE_CHECKANGLES_PLOTFIT
% Use this script to check what fit values you should restrict when reading
% in Hough data for RunEBSDrefine
% Fill in INPUT PARAMETERS section with desired parameters
% 2/23/20 (Edward Pang, MIT)

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = 'testdata_EBSDrefine/';  % path within EMdatapathname to find .ang file
inputfile = 'Scan3_Mod.ang';     % .ang file name of Hough data (located in path)
phaseid = 2;    % which phase to consider (as numbered in header of ang file)

% refined grid settings for interpolation
fitmax = 0.5;  % maximum fit value for that euler angle to be considered

% Map parameters
plotoption = 5;     % what variable to color map with (1=phi1, 2=PHI, 3=phi2, 4=CI, 5=IQ, 6=fit)
markersize = 3;     % use smaller value if you have more data points (3 is a reasonable value for 20000 map points)



%%% Parameters you don't need to change often %%%
% Paths for this computer
homepath = '/home/jonathan/epang/EMsoftfiles/EMdata/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %



% Read in .ang file
fullpath = fullfile(homepath,path,inputfile);
[euler, x, y, IQ, CI, phase, fit, phaseinfo, grid, ~] = loadang(fullpath, 1);

N = length(x);  % number of map points
xstep = grid{2};
ystep = grid{3};
ncolsodd = grid{4};
ncolseven = grid{5};
nrows = grid{6};
phi1 = euler(:,1);
PHI = euler(:,2);
phi2 = euler(:,3);


% figure out indices of phase of interest
if size(phaseinfo,1) > 1    % more than 1 phase in .ang file
    phaseid_ang = phaseid;  % for .ang files with multiple phases, phases in data rows are labeled same as header
else
    phaseid_ang = 0;    % for .ang files with one phase, phases in data rows are labeled with 0
end
indexphase = phase==phaseid_ang;  % indices of data rows of phase of interest

% get subset of x,y data for this phase
x_phase = x(indexphase);
y_phase = y(indexphase);

% pick one of these quantities for plotting
if plotoption==1
    data = phi1(indexphase);
elseif plotoption==2
    data = PHI(indexphase);
elseif plotoption==3
	data = phi2(indexphase);
elseif plotoption==4
    data = CI(indexphase);
elseif plotoption==5
    data = IQ(indexphase);
elseif plotoption==6
    data = fit(indexphase);
else
    error('Invalid plotoption input.');
end


% plot
Ncolors = 256;  % number of colors to discretize into (higher=prettier image, slower plotting)
colors = parula(Ncolors);   % rgb of colors in each row
datamin = min(data);    % min of data for color scaling
datarange = max(data)-min(data);    % range of data for color scaling


% create figure window
figure('Position',[100 100 800 round((max(y)/max(x))*800)]);
h = axes;

% loop through each color, figure out corresponding data, add to plot
for ii=1:Ncolors
    % get bounding values of data for this color
    p_low = (ii-1)/Ncolors;     % lower percentile of data to plot this iter
    p_high = ii/Ncolors;        % higher percentile of data to plot this iter
    data_low = p_low*datarange + datamin;  % min value of data to plot this iter
    data_high = p_high*datarange + datamin;  % max value of data to plot this iter
    
    % extract data for this color
    index = (data>=data_low) & (data<data_high);     % indices of points to plot
    
    % add to plot
    plot(x_phase(index),y_phase(index),'o','MarkerEdgeColor',colors(ii,:),...
        'MarkerFaceColor',colors(ii,:),'MarkerSize',markersize);
    hold on
end

% add max point
[~,imax] = max(data);
plot(x_phase(imax),y_phase(imax),'o','MarkerEdgeColor',colors(end,:),...
    'MarkerFaceColor',colors(end,:),'MarkerSize',markersize);

% tweak plot
set(gca,'Ydir','reverse');  % OIM map convention y pointing down
set(h, 'Color', 'None', 'Xtick', [], 'Ytick', []);
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);


% plot points that are ignored by the set fitmax
fit_phase = fit(indexphase);    % extract fit values for points of phase of interest
index_ignore = fit_phase>=fitmax;  % indices of data rows (in data subset of phase of interest)
plot(x_phase(index_ignore),y_phase(index_ignore),'o','MarkerEdgeColor','r',...
    'MarkerFaceColor','r','MarkerSize',markersize);     % add ignored points to plot in red


