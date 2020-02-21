%RUNPCFIT
% Use this script to perform detector fitting
% Fill in INPUT PARAMETERS section with desired parameters
% Last updated: 7/24/19 (Edward Pang, MIT)

clear

% Add paths (change these if you put snobfit and minq5 in a different location)
addpath('snobfit_v2.1/v2.1')
addpath('minq5')


%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define paths
data.path = 'CuonCu/191217_Site8/PCfit/test/';  % path within EMdatapathname to find experimental patterns and to save data

% Input data
data.img_exp_name = 'OIM Map 1_01991.png';     % name of image file located within path
data.binning = 2;        % binning mode (1, 2, 4, 8). Simulated image size will be numsx/binning x numsy/binning

% Output data
data.h5outputlogical = 1; % output h5 file? 1=yes, 0=no
data.h5output = 'test3.h5';   % Name of output file (located in path)
data.logoutputlogical = 1; % output log file of screen? 1=yes, 0=no
data.logoutput = 'test3.txt';   % Name of output file (located in path)

% Starting guesses for PC and orientation
data.L = 2508.08;    % in microns
data.xpc = 8.2329;         % in px
data.ypc = 49.1829;       % in px
data.phi1 = 112.26;      % in deg in EMsoft convention
data.PHI = 34.99;
data.phi2 = 115.50;

% Define trust region to optimize within
data.trust_L = 75;     % in microns (75 is a reasonable starting value)
data.trust_xpc = 15;     % in px   (15 is a reasonable starting value)
data.trust_ypc = 15;     % in px   (15 is a reasonable starting value)
data.trust_angle = 3;   % in degrees    (3 is a reasonable starting value)

%%% EMsoft parameters for simulating patterns and computing dot products
data.energymin = 15.0;   % energy range in the intensity summation (keV)
data.energymax = 30.0;
data.masterfile = '200106_Cu_fcc_30kV_60deg_master/Cu_fcc-master-30kV.h5';   % master pattern input file, path relative to EMdatapathname
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you know what you're doing) %


%%% Some extra parameters %%%
% detector parameters (Harvard FIB crop 87%)
data.thetac = 8;   % tilt angle of the camera (positive below horizontal, degrees)
data.delta = 7.4*480/416;   % CCD pixel size on the scintillator surface (microns)
data.numsx = 416;    % number of CCD pixels along x and y
data.numsy = 416;
data.omega = 0;      % angle between normal of sample and detector

% Paths for this computer (don't need to touch after initial setup of EMsoft)
data.homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname
data.tmppath = 'tmp/';   % path within EMdatapathname for creation of temp files

%%% Advanced parameters %%%
data.imgprebinned = 1;    % Is image already binned (matches binning below)? 1=yes, 0=no
data.h5outputlogical = 1; % output h5 file? 1=yes, 0=no
data.logoutputlogical = 1; % output log file of screen? 1=yes, 0=no
data.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
data.gammavalue = 0.33;  % gamma correction factor
data.hipassw = 0.05;    % hi pass filter w param; 0.05 is reasonable
data.nregions = 10;     % # of regions for adaptive histogram equalization
data.nthreads = 1; % for EMEBSDDI
data.maskpattern = 'y';        % use circular mask? y or n
data.r = min(data.numsx,data.numsy)/(2*data.binning);        % radius of circular mask (after binning)

%%% Define snobfit parameters
% Characteristic length scale in each direction (for scaling)
data.Lstep = 25;       % in microns
data.xpcstep = 1.5;    % in px
data.ypcstep = 1.5;    % in px
data.anglestep = 0.34;  % in degrees

% Define minimum step size for each variable
data.Lstep_min = 0.03;       % in microns
data.xpcstep_min = 0.002;   % in px
data.ypcstep_min = 0.002;   % in px
data.anglestep_min = 0.0004; % in degrees

% snobfit optimization parameters
data.ncall = 1000;   % limit on the number of function calls
data.iterstop = 100;   % stop after this number of consecutive iterations without improvement in fbest
data.npoint = 12;   % number of random start points to be generated (default=dim+6)
data.nreq = 12;     % no. of points to be generated in each call to SNOBFIT (default=dim+6)
data.p = 0.5;        % probability of generating a point of class 4 
                % (probe random region of space, global instead of local)
data.display = 'iter';  % how much info to print to screen (notify, final, off, iter)[same as fminsearch]
%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Run PC fitting
PCfit_snobfit(data)

