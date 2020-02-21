%RUNPCFIT_LOOP
% Use this script to perform detector fitting
% Fill in INPUT PARAMETERS section with desired parameters
% This version can loop through multiple images
% Original: 7/24/19 (Edward Pang, MIT)
% Change log:
% -1/8/20 ELP: Read in data and print averages, check for existing files
% initially
% -1/30/20 ELP: Instead of requiring a separate euler.txt file, directly
% input euler angles here. Save a .txt file containing the average pattern
% center from these runs.

clear

% Add paths (change these if you put snobfit and minq5 in a different location)
addpath('snobfit_v2.1/v2.1')
addpath('minq5')


%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data
data.path = 'CuonCu/191217_Site8/PCfit/test/';  % path within EMdatapathname to find experimental patterns and to save data
imagenames = [
    "OIM Map 1_01991.png"
%     "OIM Map 1_04986.png"
    ];  % name of image files to loop through (located in path)
data.binning = 2;        % binning mode (1, 2, 4, 8). Simulated image size will be numsx/binning x numsy/binning

% Define euler angles (each row phi1, PHI, phi2 in degrees)
euler = [
    112.26 34.99 115.50
%     7.88 98.64 217.81
    ];

% Output data (outputs separate file for each image)
outputprefix = 'PCfit_test3';    % output file will be [outputprefix]_[imagename(ii) without extension].h5/txt

%%% Starting guesses for PC (usually from Hough data)
data.L = 2508.08;    % in microns
data.xpc = 8.2329;         % in px
data.ypc = 49.1829;       % in px

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


N = length(imagenames);     % number of patterns to test


% check if filenames already exist
if data.h5outputlogical == 1
    for ii=1:N
        index = strfind(imagenames(ii),'.');
        s = char(imagenames(ii));
        s = s(1:index-1);
        h5output = strcat(outputprefix,'_',s,'.h5');
        h5path = fullfile(data.homepath,data.path,h5output);
        if exist(h5path,'file')==2
            error('File already exists: %s',h5path);
        end
    end
end

if data.logoutputlogical == 1
    for ii=1:N
        index = strfind(imagenames(ii),'.');
        s = char(imagenames(ii));
        s = s(1:index-1);
        logoutput = strcat(outputprefix,'_',s,'.txt');
        logpath = fullfile(data.homepath,data.path,logoutput);
        if exist(logpath,'file')==2
            error('File already exists: %s',logpath);
        end
    end
end

summarypath = strcat(data.homepath,data.path,outputprefix,'_summary.txt');
if exist(summarypath,'file')==2
    error('File already exists: %s',summarypath);
end



% loop through each image
for ii=1:N
    fprintf('\n\n=========================================================\n');
    fprintf('%s    (%g of %g)\n',imagenames(ii),ii,N);
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


% read in data
L = zeros(1,N);     % initialize
xpc = zeros(1,N);
ypc = zeros(1,N);
if data.h5outputlogical == 1
    for ii=1:N
        h5output = strcat(outputprefix,'_',s,'.h5');
        h5path = fullfile(data.homepath,data.path,h5output);
        L(ii) = h5read(h5path,'/Data/L');
        xpc(ii) = h5read(h5path,'/Data/xpc');
        ypc(ii) = h5read(h5path,'/Data/ypc');
    end
end


% compute average PC
L_ave = sum(L)/N;
L_ave_frac = L_ave/(data.numsx*data.delta);
xpc_ave = sum(xpc)/N;
xpc_ave_frac = xpc_ave/data.numsx + 0.5;
ypc_ave = sum(ypc)/N;
ypc_ave_frac = ypc_ave/data.numsy + 0.5;


% print
fprintf('\n\nFitted pattern center (average of %g patterns):\n',N);
fprintf('  xpc = %.4f (X* = %.6f)\n',xpc_ave,xpc_ave_frac);
fprintf('  ypc = %.4f (Y* = %.6f)\n',ypc_ave,ypc_ave_frac);
fprintf('  L = %.2f (Z* = %.6f)\n\n',L_ave,L_ave_frac);


% write to file
fid = fopen(summarypath,'w');
fprintf(fid,'L:');
for ii=1:N
    fprintf(fid,' %.2f,',L(ii));
end
fprintf(fid,'\nxpc:');
for ii=1:N
    fprintf(fid,' %.4f,',xpc(ii));
end
fprintf(fid,'\nypc:');
for ii=1:N
    fprintf(fid,' %.4f,',ypc(ii));
end
fprintf(fid,'\n\nFitted pattern center (average of %g patterns):\n',N);
fprintf(fid,'  xpc = %.4f (X* = %.6f)\n',xpc_ave,xpc_ave_frac);
fprintf(fid,'  ypc = %.4f (Y* = %.6f)\n',ypc_ave,ypc_ave_frac);
fprintf(fid,'  L = %.2f (Z* = %.6f)',L_ave,L_ave_frac);
fclose(fid);

fprintf('Summarized data stored in txt file: %s\n',summarypath);


