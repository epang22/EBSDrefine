function EBSDrefine_checkangles(data)
%EBSDREFINE_CHECKANGLES
% Compiles list of unique euler angles from .ang data
% Works for pseudosymmetric and non-pseudosymmetric materials
% 2/20/20 (Edward Pang, MIT)
%
%%% Input: 'data', a struct containing the following fields:
% data.homepath = '/home/jonathan/epang/EMsoftfiles/EMdata/';  % path to EMdatapathname
% data.path = '180627ELP33scan3/fullmap/';  % path within EMdatapathname to find .ang file
% data.inputfile = 'Scan3_Mod.ang';     % .ang file name of Hough data (located in path)
% 
% data.phaseid = 1;    % which phase to consider (as numbered in ang file)
% 
% % refined grid settings for interpolation
% data.angletol = 1.2;   % find all unique euler angles at least this many degrees apart as grid centers
% data.eulerminappear = 5;     % minimum number of times euler angle must appear in Hough data to be considered
% data.fitmax = 1;  % maximum fit value for that euler angle to be considered
% data.reduceeuler = 1;     % eliminate symmetry related euler angles? 1=yes, 0=no
%     % =1 is slower to find unique euler angles, but saves time later for indexing
%     % Set to 1 unless the finding euler angles step takes a really long time (>30min)
% 
% % Define pseudosymmetry variants to check (comment out if you don't want to check pseudosym)
% data.pseudosym = [
%     1 1 0 90;
%     1 -1 0 90
%     ];     % [axis_x axis_y axis_z angle(deg)]



%%% Extract input parameters from struct 'data'
homepath = data.homepath;
path = data.path;
inputfile = data.inputfile;
phaseid = data.phaseid;
angletol = data.angletol;
eulerminappear = data.eulerminappear;
fitmax = data.fitmax;
reduceeuler = data.reduceeuler;
if isfield(data,'pseudosym')
    pseudosym = data.pseudosym;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



totaltime = tic;



%%% Read in data
% read in ang file data
inputpath = fullfile(homepath,path,inputfile);   % full path of .ang file to read in
[euler,~,~,~,~,phase,fit,phaseinfo] = loadang(inputpath,1);

% get OIM symmetry ID of phase of interest
Nphases = size(phaseinfo,1);    % number of phases in .ang file
for ii=1:Nphases
    if phaseinfo{ii,1}==phaseid
        phaseindex = ii;
        symid = phaseinfo{ii,3};   % get OIM symmetry ID
        break   % found it, don't need to keep checking
    end
end

% figure out how phase of interest is labeled in ang data
if Nphases>1    % if multiple phases, ang labels them 1, 2, etc.
    phaseid_ang = phaseid;
else
    phaseid_ang = 0;    % if single phase, ang labels them 0
end


%%% Compute pseudosymmetry rotation matrices, if specified
if exist('pseudosym','var')
    Npseudo = size(pseudosym,1);    % number of pseudosym vars
    Rrot = cell(1,Npseudo);   % initialize
    for kk=1:Npseudo
        Rrot{kk} = axang2rotm([pseudosym(kk,1:3) pseudosym(kk,4)*pi/180]);
    end
end



%%% Loop through Euler angles and compile list of unique Euler angles
disp('Finding unique Euler angles...');

eulerlist1 = [];    % keep track of all unique euler angles
eulertally = [];    % keep track how many times euler angle appears (within tol)

N = size(euler,1);     % number of map points

% loop through each map point
for ii=1:N
    if phase(ii)==phaseid_ang   % only for phase of interest
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
N2 = size(eulerlist1,1);    % number of angles in eulerlist1 (picked out by misorientation)


% remove symmetry related euler angles, if specified
if reduceeuler==1
    disp('Removing symmetry related Euler angles...');

    eulerlist1_disorient = [];
    eulertally_disorient = [];

    % loop through eulerlist1
    if exist('pseudosym','var')
        for ii=1:N2
            euler1 = eulerlist1(ii,1:3);   % extract euler angle for this iter
            eulertally_temp = eulertally(ii);   % extract number of map points for this euler angle

            % compare to existing angles in eulerlist1_disorient
            counter = 1;    % if counter=1 at end, then add euler angle to list
            counter2 = 1;   % if counter2=0 at any point, then stop checking eulerlist and move on to next angle      
            for jj=1:size(eulerlist1_disorient,1)
                euler2 = eulerlist1_disorient(jj,:);   % extract euler angles for this iteration
                [~,disorient] = disorientation(euler1,euler2,symid);

                if disorient < angletol
                    counter = 0;    % don't add eulerangle to list
                    eulertally_disorient(jj) = eulertally_disorient(jj) + eulertally_temp;
                    break   % already found a match, don't need to check rest of eulerlist1_disorient
                end
                
                % pseudosym variants
                for kk=1:Npseudo
                    Reuler2 = eu2om(euler2*pi/180);     % orientation matrix of eulerlist angle
                    R2p = Rrot{kk}*Reuler2;   % orientation matrix of pseudosym var of eulerlist angle
                    euler2p = om2eu(R2p)*180/pi;   % euler angles of pseudosymmetry variant
                    [~,disorient] = disorientation(euler1,euler2p,symid);

                    if disorient < angletol
                        counter = 0;    % don't add eulerangle to list
                        counter2 = 0;   % already found a match, don't need to check rest of eulerlist
                        eulertally_disorient(jj) = eulertally_disorient(jj) + eulertally_temp;
                        break       % already found a match, don't need to check rest of pseudosym vars
                    end   
                end

                % need to keep checking eulerlist?
                if counter2==0
                    break   % already found a match, don't need to check rest of eulerlist
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
        for ii=1:N2
            euler1 = eulerlist1(ii,1:3);   % extract euler angle for this iter
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
    end
else
    % No additional disorientation filtering
    %  copy variables
    eulerlist1_disorient = eulerlist1;
    eulertally_disorient = eulertally;
end



% remove angles in eulerlist that appear less than eulerminappear times
eulerlist = eulerlist1_disorient(eulertally_disorient>=eulerminappear,:);
Neuler = size(eulerlist,1);
if Neuler<1
    fprintf('Zero unique euler angles found. Reduce eulerminappear.');
end

% print some info
fprintf('\n%g unique euler angles found within angletol\n',N2);
if reduceeuler==1
    fprintf('%g unique euler angles after removing symmetry related orientations.\n',size(eulerlist1_disorient,1));
end
fprintf('%g unique euler angles appearing at least eulerminappear times.\n\n',Neuler);
    

    
% Plot
index = length(eulertally_disorient):-1:1;
eulerminappear_line = eulerminappear*ones(size(index));

figure; plot(index,sort(eulertally_disorient));
hold on
plot(index,eulerminappear_line,'--k')
xlabel('Index of unique euler angle')
ylabel('Number of map points')

if reduceeuler==1
    t = sprintf('Filtered by disorientation < %g deg',angletol);
else
    t = sprintf('Filtered by misorientation < %g deg',angletol);
end
title(t);




toc(totaltime)

