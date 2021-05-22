function [euler, x, y, IQ, CI, phase, fit, phaseinfo, grid, PC] = loadang(inputpath, iconvert)
%LOADANG
% Read in Hough data from .ang file (multiple phases) and output euler
% angles for each phase
%%% Inputs:
% -inputpath: full path of .ang file to read in
% -iconvert: convert euler angles to EMsoft convention? 1=yes, 0=no
%%% Outputs:
% -euler: Nx3 array containing phi1, PHI, phi2 (all in deg)
% -x: Nx1 array of x position
% -y: Nx1 array of y position
% -IQ: Nx1 array of IQ
% -CI: Nx1 array of CI
% -phase: Nx1 array of phase ID
% -fit: Nx1 array of fit
% -phaseinfo: mx4 cell array. each row is phase, col1=phase id, col2=materialname, col3=symmetry id, col4: lattice param [a b c alpha beta gamma] (A, deg)
% -grid: 1x6 cell array [gridtype xstep ystep ncols_odd ncols_even nrows]
% -PC: 1x4 array [xstar ystar zstar WD]
% Original: 1/16/20 (Edward Pang, MIT)
% Change log:
% -4/20/21 ELP: make compatible with new .ang headers



% Read in data
fileID = fopen(inputpath);
C = textscan(fileID,'%s','Delimiter','\r');
fclose(fileID);


%%% Loop through each line to find start of map data
nh = regexp(C{1}{1},'# HEADER: Start','once');    % =1 if file starts with '# HEADER: Start', empty otherwise
if ~isempty(nh)
    % New .ang files, where data starts line after # HEADER: End (and first line of .ang file is # HEADER: Start)
    for ii=1:length(C{1})
        if ~isempty(regexp(C{1}{ii},'# HEADER: End','once'))
            datastartindex = ii+1;  % index of C where data begins
            break
        end
    end
else
    % Old .ang files, where data starts a couple lines after # SCANID
    for ii=1:length(C{1})
        if ~isempty(regexp(C{1}{ii},'SCANID','once'))
            datastartindex = ii+2;  % index of C where data begins
            break
        end
    end
end

N = length(C{1}) - datastartindex + 1;  % number of map points
offset = datastartindex - 1;    % difference in array index between C and my stored data arrays below



% Figure out phases, their symmetry, and their lattice params
phaseinfo = {};  % each row is phase, {1}=phase id, {2}: material name, {3}=symmetry id, {4}: lattice param [a b c alpha beta gamma] (A, deg)
counter = 0;    % number of phases
for ii=1:length(C{1})
    np = regexp(C{1}{ii},'Phase','once');    % if not empty, equals starting index of 'Phase' in string
    nn = regexp(C{1}{ii},'MaterialName','once');   % if not empty, equals starting index of 'MaterialName' in string
    ns = regexp(C{1}{ii},'Symmetry','once');    % if not empty, equals starting index of 'Symmetry' in string
    nlp = regexp(C{1}{ii},'LatticeConstants','once');    % if not empty, equals starting index of 'Symmetry' in string
    
    % get phase id
    if ~isempty(np)
        counter = counter+1;
        wholestring = C{1}{ii};
        phaseinfo{counter,1} = str2double(wholestring(np+6:end));  % store phase ID
    end

    % get material name
    if ~isempty(nn)
        wholestring = C{1}{ii};
        tabindex = regexp(wholestring,'\t','once');
        phaseinfo{counter,2} = wholestring(tabindex+1:end);  % store material name
    end
    
    % get symmetry id
    if ~isempty(ns)
        wholestring = C{1}{ii};     % whole line as single string
        wholestring2 = textscan(wholestring,'%s %s %f');    % space delimited into cell array
        phaseinfo{counter,3} = wholestring2{3};
    end

    % get lattice parameters 
    if ~isempty(nlp)
        wholestring = C{1}{ii};     % whole line as single string
        wholestring2 = textscan(wholestring,'%s %s %f %f %f %f %f %f');    % space delimited into cell array
        phaseinfo{counter,4} = [wholestring2{3} wholestring2{4} wholestring2{5} wholestring2{6} wholestring2{7} wholestring2{8}];
    end

    % stop once past header region
    if ii>datastartindex
        break
    end
end



% Read in grid info
for ii=1:length(C{1})
    if ~isempty(regexp(C{1}{ii},'GRID', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %s');
        grid{1} = gridtemp{3}{1};
    end
    if ~isempty(regexp(C{1}{ii},'XSTEP', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %f');
        grid{2} = gridtemp{3};
    end
    if ~isempty(regexp(C{1}{ii},'YSTEP', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %f');
        grid{3} = gridtemp{3};
    end
    if ~isempty(regexp(C{1}{ii},'NCOLS_ODD', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %f');
        grid{4} = gridtemp{3};
    end
    if ~isempty(regexp(C{1}{ii},'NCOLS_EVEN', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %f');
        grid{5} = gridtemp{3};
    end
    if ~isempty(regexp(C{1}{ii},'NROWS', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %f');
        grid{6} = gridtemp{3};
        break
    end
    
    % stop once past header region
    if ii>datastartindex
        break
    end
end

if ~exist('grid','var')
    grid = [];  % no grid info in .ang file, output empty vector
end


% Read in PC info
for ii=1:length(C{1})
    if ~isempty(regexp(C{1}{ii},'x-star', 'once'))
        dataxstar = textscan(C{1}{ii},'%s %s %f');
        xstar = dataxstar{3};
    end
    if ~isempty(regexp(C{1}{ii},'y-star', 'once'))
        dataystar = textscan(C{1}{ii},'%s %s %f');
        ystar = dataystar{3};
    end
    if ~isempty(regexp(C{1}{ii},'z-star', 'once'))
        datazstar = textscan(C{1}{ii},'%s %s %f');
        zstar = datazstar{3};
    end
    if ~isempty(regexp(C{1}{ii},'WorkingDistance', 'once'))
        datawd = textscan(C{1}{ii},'%s %s %f');
        wd = datawd{3};
        break
    end
    
    % stop once past header region
    if ii>datastartindex
        break
    end
end

if exist('xstar','var')
    PC = [xstar ystar zstar wd];   % put in single vector for output
else
    PC = [];    % no PC info in .ang file, output empty vector
end



% initialize some arrays
euler = zeros(N,3);
x = zeros(N,1);
y = zeros(N,1);
IQ = zeros(N,1);
CI = zeros(N,1);
phase = zeros(N,1);
fit = zeros(N,1);



%%% Read in map data
% Loop through each data point
for ii=1:N
    % Extract data
    datarow = textscan(C{1}{ii+offset},'%f');   % data row delimited into cell array
    phi1 = datarow{1}(1)*180/pi;
    PHI = datarow{1}(2)*180/pi;
    phi2 = datarow{1}(3)*180/pi;
    x(ii) = datarow{1}(4);
    y(ii) = datarow{1}(5);
    IQ(ii) = datarow{1}(6);
    CI(ii) = datarow{1}(7);
    phase(ii) = datarow{1}(8);
    fit(ii) = datarow{1}(10);

    % Convert euler angles to EMsoft convention if specified
    %  phase(ii)>0 needed to avoid error for unindexed points
    if iconvert==1 && phase(ii)>0
        % find out some info needed for conversion
        if size(phaseinfo,1)>1
            % compile vector of phase ID's
            phaseids = [];
            for jj=1:size(phaseinfo,1)
                phaseids(jj) = phaseinfo{jj,1};
            end
            
            % extract symID and lp for matching phase ID
            symID = phaseinfo{phaseids==phase(ii),3};
            lp = phaseinfo{phaseids==phase(ii),4};
        else
            symID = phaseinfo{1,3};    % in this case, phaseid=0 (does not match 1 in header)
            lp = phaseinfo{1,4};
        end

        % convert
        euler(ii,1:3) = converteuler_oim2emsoft([phi1 PHI phi2],symID,lp(4:6));
    else
        euler(ii,1) = phi1;
        euler(ii,2) = PHI;
        euler(ii,3) = phi2;
    end
end


end

