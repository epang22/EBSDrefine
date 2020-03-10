function makehdf5_3d(h5path, data)
%MAKEHDF5_3D
% Output data from ebsd_zro2pseudosym_onephase3d.m to hdf5 file
% 2/4/20 (Edward Pang, MIT)
%
%%% Inputs:
% -h5path: path to save hdf5 file
% -data: struct containing all info



% Input data
fid = H5F.create(h5path);
type_id = H5T.copy('H5T_C_S1');
H5T.set_size(type_id,'H5T_VARIABLE');
space_id = H5S.create_simple(1,1,1);
plist = 'H5P_DEFAULT';
gid = H5G.create(fid,'InputData',plist,plist,plist);
dset_id = H5D.create(fid,'InputData/homepath',type_id,space_id,plist);
data2 = {data.homepath};
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,data2);
dset_id = H5D.create(fid,'InputData/path',type_id,space_id,plist);
data2 = {data.path};
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,data2);
dset_id = H5D.create(fid,'InputData/inputfile',type_id,space_id,plist);
data2 = {data.inputfile};
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,data2);
dset_id = H5D.create(fid,'InputData/img_exp',type_id,space_id,plist);
data2 = {data.img_exp};
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,data2);
H5S.close(space_id); H5T.close(type_id); H5D.close(dset_id); H5G.close(gid); H5F.close(fid);

h5create(h5path,'/InputData/issquare',size(data.issquare),'Datatype','single');
h5write(h5path,'/InputData/issquare',single(data.issquare));
if isfield(data,'phaseid')
    h5create(h5path,'/InputData/phaseid',size(data.phaseid),'Datatype','single');
    h5write(h5path,'/InputData/phaseid',single(data.phaseid));
end
h5create(h5path,'/InputData/ipf_ht',size(data.ipf_ht),'Datatype','int32');
h5write(h5path,'/InputData/ipf_ht',int32(data.ipf_ht));
h5create(h5path,'/InputData/ipf_wd',size(data.ipf_wd),'Datatype','int32');
h5write(h5path,'/InputData/ipf_wd',int32(data.ipf_wd));
if isfield(data,'pseudosym')
    h5create(h5path,'/InputData/pseudosym',size(data.pseudosym),'Datatype','single');
    h5write(h5path,'/InputData/pseudosym',single(data.pseudosym));
end


% Pattern center
h5create(h5path,'/PatternCenter/L',size(data.L),'Datatype','single');
h5write(h5path,'/PatternCenter/L',single(data.L));
h5create(h5path,'/PatternCenter/xpc',size(data.xpc),'Datatype','single');
h5write(h5path,'/PatternCenter/xpc',single(data.xpc));
h5create(h5path,'/PatternCenter/ypc',size(data.ypc),'Datatype','single');
h5write(h5path,'/PatternCenter/ypc',single(data.ypc));


% Parameters   
h5create(h5path,'/Parameters/angletol',size(data.angletol),'Datatype','single');
h5write(h5path,'/Parameters/angletol',single(data.angletol));
h5create(h5path,'/Parameters/eulerminappear',size(data.eulerminappear),'Datatype','int32');
h5write(h5path,'/Parameters/eulerminappear',int32(data.eulerminappear));
h5create(h5path,'/Parameters/fitmax',size(data.fitmax),'Datatype','single');
h5write(h5path,'/Parameters/fitmax',single(data.fitmax));
h5create(h5path,'/Parameters/reduceeuler',size(data.reduceeuler),'Datatype','int32');
h5write(h5path,'/Parameters/reduceeuler',int32(data.reduceeuler));
if isfield(data,'refine')
    h5create(h5path,'/Parameters/refine',size(data.refine),'Datatype','int32');
    h5write(h5path,'/Parameters/refine',single(data.refine));
end
if isfield(data,'delta_angles_ref')
    h5create(h5path,'/Parameters/delta_angles_ref',size(data.delta_angles_ref),'Datatype','single');
    h5write(h5path,'/Parameters/delta_angles_ref',single(data.delta_angles_ref));
end
if isfield(data,'N_angles_ref')
    h5create(h5path,'/Parameters/N_angles_ref',size(data.N_angles_ref),'Datatype','int32');
    h5write(h5path,'/Parameters/N_angles_ref',int32(data.N_angles_ref));
end
h5create(h5path,'/Parameters/delta_angles',size(data.delta_angles),'Datatype','single');
h5write(h5path,'/Parameters/delta_angles',single(data.delta_angles));
h5create(h5path,'/Parameters/N_angles',size(data.N_angles),'Datatype','int32');
h5write(h5path,'/Parameters/N_angles',int32(data.N_angles));


% EMsoft parameters
h5create(h5path,'/EMsoftParameters/thetac',size(data.thetac),'Datatype','single');
h5write(h5path,'/EMsoftParameters/thetac',single(data.thetac));
h5create(h5path,'/EMsoftParameters/delta',size(data.delta),'Datatype','single');
h5write(h5path,'/EMsoftParameters/delta',single(data.delta));
h5create(h5path,'/EMsoftParameters/numsx',size(data.numsx),'Datatype','int32');
h5write(h5path,'/EMsoftParameters/numsx',int32(data.numsx));
h5create(h5path,'/EMsoftParameters/numsy',size(data.numsy),'Datatype','int32');
h5write(h5path,'/EMsoftParameters/numsy',int32(data.numsy));
h5create(h5path,'/EMsoftParameters/omega',size(data.omega),'Datatype','single');
h5write(h5path,'/EMsoftParameters/omega',single(data.omega));
h5create(h5path,'/EMsoftParameters/energymin',size(data.energymin),'Datatype','single');
h5write(h5path,'/EMsoftParameters/energymin',single(data.energymin));
h5create(h5path,'/EMsoftParameters/energymax',size(data.energymax),'Datatype','single');
h5write(h5path,'/EMsoftParameters/energymax',single(data.energymax));
h5create(h5path,'/EMsoftParameters/gammavalue',size(data.gammavalue),'Datatype','single');
h5write(h5path,'/EMsoftParameters/gammavalue',single(data.gammavalue));
h5create(h5path,'/EMsoftParameters/hipassw',size(data.hipassw),'Datatype','single');
h5write(h5path,'/EMsoftParameters/hipassw',single(data.hipassw));
h5create(h5path,'/EMsoftParameters/nregions',size(data.nregions),'Datatype','int32');
h5write(h5path,'/EMsoftParameters/nregions',int32(data.nregions));
h5create(h5path,'/EMsoftParameters/binning',size(data.binning),'Datatype','int32');
h5write(h5path,'/EMsoftParameters/binning',int32(data.binning));
h5create(h5path,'/EMsoftParameters/r',size(data.r),'Datatype','int32');
h5write(h5path,'/EMsoftParameters/r',int32(data.r));
h5create(h5path,'/EMsoftParameters/nthreads',size(data.nthreads),'Datatype','int32');
h5write(h5path,'/EMsoftParameters/nthreads',int32(data.nthreads));
h5create(h5path,'/EMsoftParameters/numdictsingle',size(data.numdictsingle),'Datatype','int32');
h5write(h5path,'/EMsoftParameters/numdictsingle',int32(data.numdictsingle));
h5create(h5path,'/EMsoftParameters/numexptsingle',size(data.numexptsingle),'Datatype','int32');
h5write(h5path,'/EMsoftParameters/numexptsingle',int32(data.numexptsingle));

fileattrib(h5path,'+w');
type_id = H5T.copy('H5T_C_S1');
H5T.set_size(type_id,'H5T_VARIABLE');
space_id = H5S.create_simple(1,1,1);
plist = 'H5P_DEFAULT';
fid = H5F.open(h5path,'H5F_ACC_RDWR',plist);
dset_id = H5D.create(fid,'EMsoftParameters/masterfile',type_id,space_id,plist);
data2 = {data.masterfile};
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,data2);
dset_id = H5D.create(fid,'EMsoftParameters/scalingmode',type_id,space_id,plist);
data2 = {data.scalingmode};
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,data2);
dset_id = H5D.create(fid,'EMsoftParameters/maskpattern',type_id,space_id,plist);
data2 = {data.maskpattern};
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,data2);
H5S.close(space_id); H5T.close(type_id); H5D.close(dset_id); H5F.close(fid);


% All data
h5create(h5path,'/AllData/Hougheuler',size(data.euler),'Datatype','single');
h5write(h5path,'/AllData/Hougheuler',single(data.euler));
h5create(h5path,'/AllData/eulerlist',size(data.eulerlist),'Datatype','single');
h5write(h5path,'/AllData/eulerlist',single(data.eulerlist));
h5create(h5path,'/AllData/eulerlist_all',size(data.eulerlist_all),'Datatype','single');
h5write(h5path,'/AllData/eulerlist_all',single(data.eulerlist_all));
if isfield(data,'eulerlist_pseudo')
    h5create(h5path,'/AllData/eulerlist_pseudo',size(data.eulerlist_pseudo),'Datatype','single');
    h5write(h5path,'/AllData/eulerlist_pseudo',single(data.eulerlist_pseudo));
end
h5create(h5path,'/AllData/eulerlist_grids_all',size(data.eulerlist_grids_all),'Datatype','single');
h5write(h5path,'/AllData/eulerlist_grids_all',single(data.eulerlist_grids_all));   
h5create(h5path,'/AllData/ieuler',size(data.ieuler),'Datatype','int32');
h5write(h5path,'/AllData/ieuler',int32(data.ieuler));
h5create(h5path,'/AllData/test_angles_1Drodrigues',size(data.test_angles_1Drodrigues),'Datatype','double');
h5write(h5path,'/AllData/test_angles_1Drodrigues',double(data.test_angles_1Drodrigues));
h5create(h5path,'/AllData/dp_all',size(data.dp_all),'Datatype','single');
h5write(h5path,'/AllData/dp_all',single(data.dp_all));
if isfield(data,'warnings') && ~isempty(data.warnings)
    h5create(h5path,'/AllData/warnings',size(data.warnings),'Datatype','int32');
    h5write(h5path,'/AllData/warnings',int32(data.warnings));
end
h5create(h5path,'/AllData/executiontime',size(data.executiontime),'Datatype','double');
h5write(h5path,'/AllData/executiontime',double(data.executiontime));
h5create(h5path,'/AllData/Maxinterppoint/euler',size(data.maxinterppoint_all(1:3,:,:)),'Datatype','single');
h5write(h5path,'/AllData/Maxinterppoint/euler',single(data.maxinterppoint_all(1:3,:,:)));
h5create(h5path,'/AllData/Maxinterppoint/dp',size(data.maxinterppoint_all(4,:,:)),'Datatype','single');
h5write(h5path,'/AllData/Maxinterppoint/dp',single(data.maxinterppoint_all(4,:,:)));
h5create(h5path,'/AllData/misorient_hough',size(data.misorient_hough),'Datatype','single');
h5write(h5path,'/AllData/misorient_hough',single(data.misorient_hough));
h5create(h5path,'/AllData/disorient_hough',size(data.disorient_hough),'Datatype','single');
h5write(h5path,'/AllData/disorient_hough',single(data.disorient_hough));


% Scan 1
h5create(h5path,'/Scan 1/EBSD/Data/Phi1',size(data.eulerbest(1,:)),'Datatype','single');
h5write(h5path,'/Scan 1/EBSD/Data/Phi1',single(data.eulerbest(1,:)));
h5create(h5path,'/Scan 1/EBSD/Data/Phi',size(data.eulerbest(2,:)),'Datatype','single');
h5write(h5path,'/Scan 1/EBSD/Data/Phi',single(data.eulerbest(2,:)));
h5create(h5path,'/Scan 1/EBSD/Data/Phi2',size(data.eulerbest(3,:)),'Datatype','single');
h5write(h5path,'/Scan 1/EBSD/Data/Phi2',single(data.eulerbest(3,:)));
h5create(h5path,'/Scan 1/EBSD/Data/CI',size(data.CI),'Datatype','single');
h5write(h5path,'/Scan 1/EBSD/Data/CI',single(data.CI));
h5create(h5path,'/Scan 1/EBSD/Data/Fit',size(data.fit),'Datatype','single');
h5write(h5path,'/Scan 1/EBSD/Data/Fit',single(data.fit));
h5create(h5path,'/Scan 1/EBSD/Data/IQ',size(data.IQ),'Datatype','single');
h5write(h5path,'/Scan 1/EBSD/Data/IQ',single(data.IQ));
h5create(h5path,'/Scan 1/EBSD/Data/Phase',size(data.phase),'Datatype','int32');
h5write(h5path,'/Scan 1/EBSD/Data/Phase',int32(data.phase));
h5create(h5path,'/Scan 1/EBSD/Data/TopDotProductList',size(data.dpbest),'Datatype','single');
h5write(h5path,'/Scan 1/EBSD/Data/TopDotProductList',single(data.dpbest));
h5create(h5path,'/Scan 1/EBSD/Data/X Position',size(data.x),'Datatype','single');
h5write(h5path,'/Scan 1/EBSD/Data/X Position',single(data.x));
h5create(h5path,'/Scan 1/EBSD/Data/Y Position',size(data.y),'Datatype','single');
h5write(h5path,'/Scan 1/EBSD/Data/Y Position',single(data.y));



