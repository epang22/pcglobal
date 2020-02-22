%RUNPCFIT
% Use this script to fit the pattern center to a single pattern
% Fill in INPUT PARAMETERS section with desired parameters
% 2/21/20 (Edward Pang, MIT)

clear

% Add paths (change these if you put snobfit and minq5 in a different location)
addpath('snobfit_v2.1/v2.1')
addpath('minq5')


%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define paths
data.path = 'testdata/';  % path within EMdatapathname to find experimental patterns and to save data

% Input data
data.img_exp_name = 'EDAX-Ni_Scan1_x126y0.png';     % name of image file located within path
data.binning = 8;        % binning mode (1, 2, 4, 8). Simulated image size will be numsx/binning x numsy/binning

% Output data
data.h5output = 'output.h5';   % Name of output file (located in path)
data.logoutput = 'output.txt';   % Name of output file (located in path)

% Starting guesses for PC and orientation
data.L = 15870.0234;    % in microns
data.xpc = 3.48576;         % in px
data.ypc = 114.20352;       % in px
data.phi1 = 35.2655;      % in deg in EMsoft convention
data.PHI = 109.96;
data.phi2 = 205.439;

% Define trust region to optimize within
data.trust_L = 500;     % in microns
data.trust_xpc = 5;     % in px
data.trust_ypc = 5;     % in px
data.trust_angle = 1;   % in degrees

% EMsoft parameters for the master pattern
data.energymin = 10.0;   % energy range in the intensity summation (keV)
data.energymax = 20.0;
data.masterfile = 'DItutorial_Ni-master-20kV/Ni-master-20kV.h5';   % master pattern input file, path relative to EMdatapathname



%%% Parameters you don't need to change often %%%
% Paths for this computer
data.homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)
data.tmppath = 'tmp/';   % path within EMdatapathname for creation of temp files (you need to manually create this directory)

% Is image already binned (matches binning below)? 1=yes, 0=no
%   Allows you to further bin down patterns you read in
data.imgprebinned = 1; 	% 1=yes (already binned, image size matches numsx/binning by numsy/binning), 0=no (you want to additionally bin)

% Files to output
data.h5outputlogical = 1; % output h5 file? 1=yes, 0=no
data.logoutputlogical = 1; % output log file of screen? 1=yes, 0=no

% Detector parameters
data.thetac = 10;   % tilt angle of the camera (positive below horizontal, degrees)
data.delta = 59.2;   % CCD pixel size on the scintillator surface (microns)
data.numsx = 480;    % number of CCD pixels along x and y
data.numsy = 480;
data.omega = 0;      % angle between normal of sample and detector

% Some EMsoft parameters for dot product computations
data.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
data.gammavalue = 0.33;  % gamma correction factor
data.hipassw = 0.05;    % hi pass filter w param; 0.05 is reasonable
data.nregions = 10;     % # of regions for adaptive histogram equalization
data.nthreads = 1; % for EMEBSDDI
data.maskpattern = 'y';        % use circular mask? y or n
data.r = min(data.numsx,data.numsy)/(2*data.binning);        % radius of circular mask (after binning)

%%% Some more SNOBFIT parameters
%   Default values are fine in most cases, can use the FindParameters.m program to find the optimal parameters for your system
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

% Optimization parameters (default values are generally fine)
data.ncall = 1000;   % limit on the number of function calls
data.iterstop = 100;   % stop after this number of consecutive iterations without improvement in fbest
data.npoint = 12;   % number of random start points to be generated (default=dim+6)
data.nreq = 12;     % no. of points to be generated in each call to SNOBFIT (default=dim+6)
data.p = 0.5;        % probability of generating a point of class 4 
                % (probe random region of space, global instead of local)
data.display = 'iter';  % how much info to print to screen (notify, final, off, iter)[same as fminsearch]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %



%%% Run PC fitting
PCfit_snobfit(data)

