%RUNPCFIT_LOOP
% Use this script to perform detector fitting
% Fill in INPUT PARAMETERS section with desired parameters
% This version can loop through multiple images
% Requires an euler.txt file (in EMsoft format) that specifies starting
% orientations of each pattern
% Last updated: 7/24/19 (Edward Pang, MIT)

clear

% Add paths (change these if you put snobfit and minq5 in a different location)
addpath('snobfit_v2.1/v2.1')
addpath('minq5')


%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define paths
data.homepath = '/home/jonathan/epang/EMsoftfiles/EMdata/';  % path to EMdatapathname
  %^This is EMsoft's EMdatapathname for your installation. If you are not sure where this is, you can find it at: ~/.config/EMsoft/EMsoftConfig.json
data.path = 'testdata/';  % path within EMdatapathname to find experimental patterns and to save data
data.tmppath = 'tmp/';   % path within EMdatapathname for creation of temp files

% Input data
imagenames = [
    "EDAX-Ni_Scan1_x126y0.png"
    "EDAX-Ni_Scan1_x126y2.png"
    "EDAX-Ni_Scan1_x126y3.png"
    ];  % name of image files to loop through (located in path)
orifile = 'euler.txt';  % file specifying starting orientations, in same order as imagenames (located in path)
data.imgprebinned = 1;    % Is image already binned (matches binning below)? 1=yes, 0=no

% Output data (outputs separate file for each image)
data.h5outputlogical = 1; % output h5 file? 1=yes, 0=no
data.logoutputlogical = 1; % output log file of screen? 1=yes, 0=no
outputprefix = 'output';    % output file will be [outputprefix]_[imagename(ii) without extension].h5/txt


%%% Starting guesses for PC
data.L = 15870.0234;    % in microns
data.xpc = 3.48576;         % in px
data.ypc = 114.20352;       % in px


%%% Define snobfit parameters
% Characteristic length scale in each direction (for scaling)
data.Lstep = 225;       % in microns
data.xpcstep = 2.25;    % in px
data.ypcstep = 2.25;    % in px
data.anglestep = 0.43;  % in degrees

% Define minimum step size for each variable
data.Lstep_min = 0.25;       % in microns
data.xpcstep_min = 0.0025;   % in px
data.ypcstep_min = 0.0025;   % in px
data.anglestep_min = 0.0005; % in degrees

% Define trust region to optimize within
data.trust_L = 500;     % in microns
data.trust_xpc = 5;     % in px
data.trust_ypc = 5;     % in px
data.trust_angle = 1;   % in degrees

% snobfit optimization parameters
data.ncall = 20;   % limit on the number of function calls
data.iterstop = 100;   % stop after this number of consecutive iterations without improvement in fbest
data.npoint = 12;   % number of random start points to be generated (default=dim+6)
data.nreq = 12;     % no. of points to be generated in each call to SNOBFIT (default=dim+6)
data.p = 0.5;        % probability of generating a point of class 4 
                % (probe random region of space, global instead of local)
data.display = 'iter';  % how much info to print to screen (notify, final, off, iter)[same as fminsearch]

%%% EMsoft parameters for simulating patterns and computing dot products
data.thetac = 10;   % tilt angle of the camera (positive below horizontal, degrees)
data.delta = 59.2;   % CCD pixel size on the scintillator surface (microns)
data.numsx = 480;    % number of CCD pixels along x and y
data.numsy = 480;
data.omega = 0;      % angle between normal of sample and detector
data.energymin = 10.0;   % energy range in the intensity summation (keV)
data.energymax = 20.0;
data.masterfile = 'DItutorial_Ni-master-20kV/Ni-master-20kV.h5';   % master pattern input file, path relative to EMdatapathname
data.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
data.gammavalue = 0.33;  % gamma correction factor
data.maskpattern = 'y';        % use circular mask? y or n
data.hipassw = 0.05;    % hi pass filter w param; 0.05 is reasonable
data.nregions = 10;     % # of regions for adaptive histogram equalization
data.binning = 8;        % binning mode (1, 2, 4, 8). Simulated image size will be numsx/binning x numsy/binning
data.r = min(data.numsx,data.numsy)/(2*data.binning);        % radius of circular mask (after binning)
data.nthreads = 8; % for EMEBSDDI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %


% load orientation file
euler = dlmread(fullfile(data.homepath,data.path,orifile),' ',2,1);


% loop through each image
for ii=1:length(imagenames)
    
    fprintf('\n\n=========================================================\n');
    fprintf('%s    (%g of %g)\n',imagenames(ii),ii,length(imagenames));
    fprintf('=========================================================\n');

    data.img_exp_name = imagenames(ii);
    data.phi1 = euler(ii,1);
    data.PHI = euler(ii,2);
    data.phi2 = euler(ii,3);

    % figure out names of output files
    index = strfind(imagenames(ii),'.');
    s = char(imagenames(ii));
    s = s(1:index-1);
    data.h5output = strcat(outputprefix,'_',s,'.h5');
    data.logoutput = strcat(outputprefix,'_',s,'.txt');

    
    %%% Run PC fitting
    PCfit_snobfit(data)

end

