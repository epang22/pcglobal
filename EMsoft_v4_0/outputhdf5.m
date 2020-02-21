function outputhdf5(h5path,data,PC_final)
%OUTPUTHDF5 Take data and output as hdf5 file
%%% Inputs:
% -h5path: string containing full path of .h5 file to create and save data in
% -data: struct containing parameters to output
% -PC_final: 1x7 array containing optimized PC: L, xpc, ypc, phi1, PHI, phi2, dp
%
% Last updated: 7/24/19 (Edward Pang, MIT)


    %%% Extract values from struct 'data'
    L = data.L;
    xpc = data.xpc;
    ypc = data.ypc;
    phi1 = data.phi1;
    PHI = data.PHI;
    phi2 = data.phi2;
    Lstep = data.Lstep;
    xpcstep = data.xpcstep;
    ypcstep = data.ypcstep;
    anglestep = data.anglestep;
    Lstep_min = data.Lstep_min;
    xpcstep_min = data.xpcstep_min;
    ypcstep_min = data.ypcstep_min;
    anglestep_min = data.anglestep_min;
    trust_L = data.trust_L;
    trust_xpc = data.trust_xpc;
    trust_ypc = data.trust_ypc;
    trust_angle = data.trust_angle;
    ncall = data.ncall;
    iterstop = data.iterstop;
    npoint = data.npoint;
    nreq = data.nreq;
    p = data.p;
    thetac = data.thetac;
    delta = data.delta;
    numsx = data.numsx;
    numsy = data.numsy;
    omega = data.omega;
    energymin = data.energymin;
    energymax = data.energymax;
    masterfile = data.masterfile;
    binning = data.binning;
    scalingmode = data.scalingmode;
    gammavalue = data.gammavalue;
    maskpattern = data.maskpattern;
    nthreads = data.nthreads;
    hipassw = data.hipassw;
    nregions = data.nregions;
    r = data.r;

    % Write starting guesses for PC and orientation
    h5create(h5path,'/StartingGuess/L',size(L),'Datatype','double');
    h5write(h5path,'/StartingGuess/L',L);
    h5create(h5path,'/StartingGuess/xpc',size(xpc),'Datatype','double');
    h5write(h5path,'/StartingGuess/xpc',xpc);
    h5create(h5path,'/StartingGuess/ypc',size(ypc),'Datatype','double');
    h5write(h5path,'/StartingGuess/ypc',ypc);
    h5create(h5path,'/StartingGuess/phi1',size(phi1),'Datatype','double');
    h5write(h5path,'/StartingGuess/phi1',phi1);
    h5create(h5path,'/StartingGuess/PHI',size(PHI),'Datatype','double');
    h5write(h5path,'/StartingGuess/PHI',PHI);
    h5create(h5path,'/StartingGuess/phi2',size(phi2),'Datatype','double');
    h5write(h5path,'/StartingGuess/phi2',phi2);

    % Write snobfit parameters
    h5create(h5path,'/SnobfitParameters/Lstep',size(Lstep),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/Lstep',Lstep);
    h5create(h5path,'/SnobfitParameters/xpcstep',size(xpcstep),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/xpcstep',xpcstep);
    h5create(h5path,'/SnobfitParameters/ypcstep',size(ypcstep),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/ypcstep',ypcstep);
    h5create(h5path,'/SnobfitParameters/anglestep',size(anglestep),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/anglestep',anglestep);
    h5create(h5path,'/SnobfitParameters/Lstep_min',size(Lstep_min),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/Lstep_min',Lstep_min);
    h5create(h5path,'/SnobfitParameters/xpcstep_min',size(xpcstep_min),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/xpcstep_min',xpcstep_min);
    h5create(h5path,'/SnobfitParameters/ypcstep_min',size(ypcstep_min),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/ypcstep_min',ypcstep_min);
    h5create(h5path,'/SnobfitParameters/anglestep_min',size(anglestep_min),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/anglestep_min',anglestep_min);    
    h5create(h5path,'/SnobfitParameters/trust_L',size(trust_L),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/trust_L',trust_L);
    h5create(h5path,'/SnobfitParameters/trust_xpc',size(trust_xpc),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/trust_xpc',trust_xpc);
    h5create(h5path,'/SnobfitParameters/trust_ypc',size(trust_ypc),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/trust_ypc',trust_ypc);
    h5create(h5path,'/SnobfitParameters/trust_angle',size(trust_angle),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/trust_angle',trust_angle);
    h5create(h5path,'/SnobfitParameters/ncall',size(ncall),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/ncall',ncall);     
    h5create(h5path,'/SnobfitParameters/iterstop',size(iterstop),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/iterstop',iterstop);     
    h5create(h5path,'/SnobfitParameters/npoint',size(npoint),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/npoint',npoint); 
    h5create(h5path,'/SnobfitParameters/nreq',size(nreq),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/nreq',nreq); 
    h5create(h5path,'/SnobfitParameters/p',size(p),'Datatype','double');
    h5write(h5path,'/SnobfitParameters/p',p); 

    % Write EMsoft parameters
    h5create(h5path,'/EMsoftParameters/thetac',size(thetac),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/thetac',thetac); 
    h5create(h5path,'/EMsoftParameters/delta',size(delta),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/delta',delta); 
    h5create(h5path,'/EMsoftParameters/numsx',size(numsx),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/numsx',numsx); 
    h5create(h5path,'/EMsoftParameters/numsy',size(numsy),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/numsy',numsy); 
    h5create(h5path,'/EMsoftParameters/omega',size(omega),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/omega',omega); 
    h5create(h5path,'/EMsoftParameters/energymin',size(energymin),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/energymin',energymin); 
    h5create(h5path,'/EMsoftParameters/energymax',size(energymax),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/energymax',energymax);
    
    file_id = H5F.open(h5path,'H5F_ACC_RDWR','H5P_DEFAULT');
    filetype = H5T.copy ('H5T_C_S1');
    H5T.set_size (filetype,length(masterfile));
    memtype = H5T.copy ('H5T_C_S1');
    H5T.set_size (memtype,length(masterfile));
    space_id = H5S.create_simple(1,1,1);
    dataset_id = H5D.create(file_id, '/EMsoftParameters/masterfile', filetype, space_id, 'H5P_DEFAULT');
    H5D.write (dataset_id, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', masterfile);
    
    H5T.set_size (filetype,length(scalingmode));
    H5T.set_size (memtype,length(scalingmode));
    dataset_id = H5D.create(file_id, '/EMsoftParameters/scalingmode', filetype, space_id, 'H5P_DEFAULT');
    H5D.write (dataset_id, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', scalingmode);
 
    h5create(h5path,'/EMsoftParameters/gammavalue',size(gammavalue),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/gammavalue',gammavalue); 
    
    H5T.set_size (filetype,length(maskpattern));
    H5T.set_size (memtype,length(maskpattern));
    dataset_id = H5D.create(file_id, '/EMsoftParameters/maskpattern', filetype, space_id, 'H5P_DEFAULT');
    H5D.write (dataset_id, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', maskpattern);

    h5create(h5path,'/EMsoftParameters/hipassw',size(hipassw),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/hipassw',hipassw); 
    h5create(h5path,'/EMsoftParameters/nregions',size(nregions),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/nregions',nregions); 
    h5create(h5path,'/EMsoftParameters/binning',size(binning),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/binning',binning); 
    h5create(h5path,'/EMsoftParameters/r',size(r),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/r',r); 
    h5create(h5path,'/EMsoftParameters/nthreads',size(nthreads),'Datatype','double');
    h5write(h5path,'/EMsoftParameters/nthreads',nthreads); 
    
    % Write final PC
    h5create(h5path,'/Data/L',size(PC_final(1)),'Datatype','double');
    h5write(h5path,'/Data/L',PC_final(1));
    h5create(h5path,'/Data/xpc',size(PC_final(2)),'Datatype','double');
    h5write(h5path,'/Data/xpc',PC_final(2));
    h5create(h5path,'/Data/ypc',size(PC_final(3)),'Datatype','double');
    h5write(h5path,'/Data/ypc',PC_final(3));
    h5create(h5path,'/Data/dp',size(PC_final(7)),'Datatype','double');
    h5write(h5path,'/Data/dp',PC_final(7));
    h5create(h5path,'/Data/euler',size(PC_final(4:6)),'Datatype','double');
    h5write(h5path,'/Data/euler',PC_final(4:6));
    
    % Close hdf5 file
    H5D.close(dataset_id);
    H5S.close(space_id);
    H5T.close(filetype);
    H5T.close(memtype);
    H5F.close(file_id);

end

