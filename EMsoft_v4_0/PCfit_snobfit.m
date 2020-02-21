function PCfit_snobfit( data )
%PCFIT_SNOBFIT Fit PC for a single pattern by snobfit in 6d space
%%% Input: 'data', a struct containing the following fields:
% % Define paths
% data.homepath = '/home/jonathan/epang/EMsoftfiles/EMdata/';  % path to EMdatapathname
% data.path = '190621_DItutorial_Ni/testdata/';  % path within EMdatapathname to find experimental patterns and to save data
% data.tmppath = 'tmp/';   % path within EMdatapathname for creation of temp files
% 
% % Input data
% data.img_exp_name = 'EDAX-Ni_Scan1_x126y0.png';     % name of image file located within path
% data.imgprebinned = 1;    % Is image already binned (matches binning below)? 1=yes, 0=no
% 
% % Output data
% data.h5outputlogical = 1; % output h5 file? 1=yes, 0=no
% data.h5output = 'output.h5';   % Name of output file (located in path)
% data.logoutputlogical = 1; % output log file of screen? 1=yes, 0=no
% data.logoutput = 'output.txt';   % Name of output file (located in path)
% 
% %%% Starting guesses for PC and orientation
% data.L = 15870.0234;    % in microns
% data.xpc = 3.48576;         % in px
% data.ypc = 114.20352;       % in px
% data.phi1 = 35.27;      % in deg in EMsoft convention
% data.PHI = 109.96;
% data.phi2 = 205.44;
% 
% %%% Define snobfit parameters
% % Characteristic length scale in each direction (for scaling)
% data.Lstep = 225;       % in microns
% data.xpcstep = 2.25;    % in px
% data.ypcstep = 2.25;    % in px
% data.anglestep = 0.43;  % in degrees
% 
% % Define minimum step size for each variable
% data.Lstep_min = 0.25;       % in microns
% data.xpcstep_min = 0.0025;   % in px
% data.ypcstep_min = 0.0025;   % in px
% data.anglestep_min = 0.0005; % in degrees
% 
% % Define trust region to optimize within
% data.trust_L = 500;     % in microns
% data.trust_xpc = 5;     % in px
% data.trust_ypc = 5;     % in px
% data.trust_angle = 1;   % in degrees
% 
% % snobfit optimization parameters
% data.ncall = 1000;   % limit on the number of function calls
% data.iterstop = 100;   % stop after this number of consecutive iterations without improvement in fbest
% data.npoint = 12;   % number of random start points to be generated (default=dim+6)
% data.nreq = 12;     % no. of points to be generated in each call to SNOBFIT (default=dim+6)
% data.p = 0.5;        % probability of generating a point of class 4 
%                 % (probe random region of space, global instead of local)
% data.display = 'iter';  % how much info to print to screen (notify, final, off, iter)[same as fminsearch]
% 
% %%% EMsoft parameters for simulating patterns and computing dot products
% data.thetac = 10;   % tilt angle of the camera (positive below horizontal, degrees)
% data.delta = 59.2;   % CCD pixel size on the scintillator surface (microns)
% data.numsx = 480;    % number of CCD pixels along x and y
% data.numsy = 480;
% data.omega = 0;      % angle between normal of sample and detector
% data.energymin = 10.0;   % energy range in the intensity summation (keV)
% data.energymax = 20.0;
% data.masterfile = 'DItutorial_Ni-master-20kV/Ni-master-20kV.h5';   % master pattern input file, path relative to EMdatapathname
% data.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
% data.gammavalue = 0.33;  % gamma correction factor
% data.maskpattern = 'y';        % use circular mask? y or n
% data.hipassw = 0.05;    % hi pass filter w param; 0.05 is reasonable
% data.nregions = 10;     % # of regions for adaptive histogram equalization
% data.binning = 8;        % binning mode (1, 2, 4, 8). Simulated image size will be numsx/binning x numsy/binning
% data.r = min(data.numsx,data.numsy)/(2*data.binning);        % radius of circular mask (after binning)
% data.nthreads = 8; % for EMEBSDDI

% Original: 7/24/19 (Edward Pang, MIT)
% Change log:
% -11/16/19 ELP: fix bug, rename output file if exists
% -12/25/19 ELP: change euler search edge warning criterion


%%% Extract input parameters from struct 'data'
homepath = data.homepath;
path = data.path;
tmppath = data.tmppath;
img_exp_name = data.img_exp_name;
h5outputlogical = data.h5outputlogical;
h5output = data.h5output;
logoutputlogical = data.logoutputlogical;
logoutput = data.logoutput;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


totaltime = tic;

% clear tmp folder
[~,~] = system(sprintf('rm %s%s*',homepath,tmppath));     % delete all files

% define paths
h5path = fullfile(homepath,path,h5output);    % full path to h5 output file
logpath = fullfile(homepath,path,logoutput);    % full path to txt log output file
img_exp_path = fullfile(homepath,path,img_exp_name);    % full path to pattern image file


% check to see if image exists, if so continue, if not print message
if exist(img_exp_path,'file')==2
    
    % Check if output files already exist. If so, rename
    if exist(h5path,'file')==2 || exist(logpath,'file')==2
        while exist(h5path,'file')==2 || exist(logpath,'file')==2
            h5output = [h5output(1:end-3) '_1.h5'];
            h5path = fullfile(homepath,path,h5output);    % full path to h5 output file
            logoutput = [logoutput(1:end-4) '_1.txt'];
            logpath = fullfile(homepath,path,logoutput);    % full path to log output file
        end
        fprintf('h5 output file renamed to: %s\n',h5path);
        fprintf('log output file renamed to: %s\n',logpath);
    end
    
    % open log diary file
    if logoutputlogical==1
        diary(sprintf('%s',logpath))
    end
    
    
    % load experimental pattern, bin it down, save to binary .data file
    outname = 'test.data';  % name of binary .data file
    data.img_exp = outname;     % add to data, for emebsddi_dp.m and loadpattern.m
    loadpattern(data);


    % Compute scale factors for optimization
    Lfactor = Lstep/0.00025;
    xpcfactor = xpcstep/0.00025;
    ypcfactor = ypcstep/0.00025;
    Rfactor = tand(anglestep/2)/0.00025;
    
    
    % Compute snobfit resolution vector (minimum step size in x that is meaningful)
    dx = [Lstep_min xpcstep_min ypcstep_min anglestep_min anglestep_min anglestep_min];
    dx_scale = [dx(1)/Lfactor dx(2)/xpcfactor dx(3)/ypcfactor tand(dx(4)/2)/Rfactor ...
        tand(dx(5)/2)/Rfactor tand(dx(6)/2)/Rfactor]';   % rescaled dx for snobfit

    
    % Compute bounds
    anglebound = tand(trust_angle/2)/Rfactor;
    ub = [trust_L/Lfactor trust_xpc/xpcfactor trust_ypc/ypcfactor anglebound anglebound anglebound]';    % upper bound
    lb = -ub;    % lower bound
    
    % optimize L,xpc,ypc,euler angles by calling snobfit
    q1 = eu2qu([phi1 PHI phi2]*pi/180);     % Quaternion of starting Euler angles (center of grid)
    fun = @(x)-1*emebsddi_wrapper(x,L,xpc,ypc,q1,Lfactor,xpcfactor,ypcfactor,Rfactor,data);   % define function to minimize
    [testminparam,fval] = ebsd_snobfit(fun,dx_scale,lb,ub,data);    % Snobfit

    % store results
    PC_final = zeros(1,7);    % Initialize array to store optimized PC: L, xpc, ypc, phi1, PHI, phi2, dp
    PC_final(1) = testminparam(1)*Lfactor + L;
    PC_final(2) = testminparam(2)*xpcfactor + xpc;
    PC_final(3) = testminparam(3)*ypcfactor + ypc;
    [euler, qrot] = rodmis2eu(testminparam(4:6)*Rfactor,q1);    % convert Rodrigues misorientation to euler angles
    PC_final(4) = euler(1);
    PC_final(5) = euler(2);
    PC_final(6) = euler(3);
    PC_final(7) = -fval;           

    

    % print optimized PC to screen
    fprintf('Optimized pattern center and orientation:\n');
    fprintf('  L = %.2f\n',PC_final(1));
    fprintf('  xpc = %.4f\n',PC_final(2));
    fprintf('  ypc = %.4f\n',PC_final(3));
    fprintf('  phi1 = %.2f deg\n',PC_final(4));
    fprintf('  PHI = %.2f deg\n',PC_final(5));
    fprintf('  phi2 = %.2f deg\n',PC_final(6));
    fprintf('  dot product = %.6f\n\n',PC_final(7));



    % Check if final solution is near edge of search radius, if so print warning
    if abs(PC_final(1)-L)>0.95*trust_L
        fprintf('WARNING: L is near search edge. trust_L = %.2f, delta_L = %.2f.\n',trust_L,abs(PC_final(1)-L));
    end
    if abs(PC_final(2)-xpc)>0.95*trust_xpc
        fprintf('WARNING: xpc is near search edge. trust_xpc = %.4f, delta_xpc = %.4f.\n',trust_xpc,abs(PC_final(2)-xpc));
    end
    if abs(PC_final(3)-ypc)>0.95*trust_ypc
        fprintf('WARNING: ypc is near search edge. trust_ypc = %.4f, delta_ypc = %.4f.\n',trust_ypc,abs(PC_final(3)-ypc));
    end
    % old version, not correct because it only checks corners of boxes, not edges/faces
%     if 2*acosd(dot(q1,qrot))>0.95*sqrt(3)*trust_angle
%         fprintf('WARNING: euler is near search edge. sqrt(3)*trust_angle = %.4f, delta_angle = %.4f.\n',sqrt(3)*trust_angle,2*acosd(dot(q1,qrot)));
%     end

    % new version
    if testminparam(4)>0.95*anglebound || testminparam(4)<-0.95*anglebound || testminparam(5)>0.95*anglebound || testminparam(5)<-0.95*anglebound || testminparam(6)>0.95*anglebound || testminparam(6)<-0.95*anglebound
        fprintf('WARNING: euler is near search edge.\n');
    end

    fprintf('\n');
    toc(totaltime)
    


    % save data to hdf5 file
    if h5outputlogical==1
        outputhdf5(h5path,data,PC_final);
    end


% print message to screen if image not found
else
    fprintf('\n===========================================================\n');
    fprintf('  File not found: %s\n',img_exp_path);
    fprintf('\n===========================================================\n');
end


% turn diary log off
if logoutputlogical==1
    diary off
end


end

