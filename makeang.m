function makeang(outpath, data, iconvert)
%MAKEANG
% Make .ang file from data
% 1/30/20 (Edward Pang, MIT)
%
%%% Inputs:
% -outpath: full path to export .ang file to
% -data: struct containing data to export
% -iconvert: convert euler angles to OIM convention? 1=yes, 0=no



% Extract data
eulerbest = data.eulerbest;     % Note this is 3xN in rad
CI = data.CI;
fit = data.fit;
IQ = data.IQ;
phase = data.phase;
x = data.x;
y = data.y;
phaseinfo = data.phaseinfo;
% info to print in info line
if isfield(data,'info')
    info = data.info;
else
    info = '';
end
% grid info: 1x6 cell array {gridtype xstep ystep ncols_odd ncols_even nrows}
if isfield(data,'grid')
    grid = data.grid;
end
% PC info: 1x4 array [xstar ystar zstar WD]
if isfield(data,'PC')
    PC = data.PC;
end

N = length(x);  % number of map points


%%% write to file
fileid = fopen(outpath,'w');     % Open file

% PC info, if present
if exist('PC','var')
    fprintf(fileid,'# TEM_PIXperUM:          1.000000\n');
    fprintf(fileid,'# x-star:                %.6f\n',PC(1));
    fprintf(fileid,'# y-star:                %.6f\n',PC(2));
    fprintf(fileid,'# z-star:                %.6f\n',PC(3));
    fprintf(fileid,'# WorkingDistance:       %.6f\n',PC(4));
    fprintf(fileid,'#\n');
end

% phase info
for ii=1:size(phaseinfo,1)
    fprintf(fileid,'# Phase %g\n',phaseinfo{ii,1});
    fprintf(fileid,'# MaterialName  \t%s\n',phaseinfo{ii,2});
    fprintf(fileid,'# Info  \t%s\n',info);
    fprintf(fileid,'# Symmetry              %.0f\n',phaseinfo{ii,3});
    fprintf(fileid,'# LatticeConstants      %5.3f %5.3f %5.3f %7.3f %7.3f %7.3f\n',phaseinfo{ii,4}(1),phaseinfo{ii,4}(2),phaseinfo{ii,4}(3),phaseinfo{ii,4}(4),phaseinfo{ii,4}(5),phaseinfo{ii,4}(6));
    fprintf(fileid,'# NumberFamilies        0\n');
    fprintf(fileid,'#\n');
end

% grid info, if present
if exist('grid','var')
    fprintf(fileid,'# GRID: %s\n',grid{1});
    fprintf(fileid,'# XSTEP: %.6f\n',grid{2});
    fprintf(fileid,'# YSTEP: %.6f\n',grid{3});
    fprintf(fileid,'# NCOLS_ODD: %g\n',grid{4});
    fprintf(fileid,'# NCOLS_EVEN: %g\n',grid{5});
    fprintf(fileid,'# NROWS: %g\n',grid{6});
    fprintf(fileid,'#\n');
end

fprintf(fileid,'# SCANID: \t\n');
fprintf(fileid,'#\n');


% Convert Euler angles to OIM convention if specified
if iconvert == 1
    % initialize
    phi1 = zeros(size(x));
    PHI = zeros(size(x));
    phi2 = zeros(size(x));
    
    % loop through each map point
    for ii=1:N
        % extract euler angles for this point, convert to deg
        euler = eulerbest(:,ii)'*180/pi;
        
        % find out some info needed for conversion
        if size(phaseinfo,1)>1
            % compile vector of phase ID's
            phaseids = [];
            for jj=1:size(phaseinfo,1)
                phaseids(jj) = phaseinfo{jj,1};
            end
            
            % extract symID and lp for matching phase ID
            symID = phaseinfo{phaseids==phase(ii),3};
            lpangles = phaseinfo{phaseids==phase(ii),4}(4:6);
        else
            symID = phaseinfo{1,3};
            lpangles = phaseinfo{1,4}(4:6);
        end

        % convert euler angles
        euler_oim = converteuler_emsoft2oim(euler, symID, lpangles);
        phi1(ii) = euler_oim(1)*pi/180;     % in rad
        PHI(ii) = euler_oim(2)*pi/180;
        phi2(ii) = euler_oim(3)*pi/180;
    end
else
    phi1 = eulerbest(1,:)';
    PHI = eulerbest(2,:)';
    phi2 = eulerbest(3,:)';
end


% Print data to file
for ii = 1:N
%     fprintf(fileid,'%9.5f %9.5f %9.5f %12.5f %12.5f %.1f %6.3f  %.0f      1  %.3f\n',phi1(ii),PHI(ii),phi2(ii),x(ii),y(ii),IQ(ii),CI(ii),phase(ii),fit(ii));    % OIM format
    fprintf(fileid,'%9.5f %9.5f %9.5f %12.5f %12.5f %.1f %8.5f  %.0f      1  %.5f\n',phi1(ii),PHI(ii),phi2(ii),x(ii),y(ii),IQ(ii),CI(ii),phase(ii),fit(ii));    % more decimal places
end


fclose('all');

