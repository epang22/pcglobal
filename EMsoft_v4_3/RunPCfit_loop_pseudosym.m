%RUNPCFIT_LOOP_PSEUDOSYM
% Use this script to fit pattern centers to multiple patterns for
% pseudosymmetric materials. This code will fit to all pseudosymmetry
% variants if supplied with the pseudosymmetry rotation.
% Fill in INPUT PARAMETERS section with desired parameters
% This version can loop through multiple images
% 5/8/21 (Edward Pang, MIT)

clear

% Add paths (change these if you put snobfit and minq5 in a different location)
addpath('snobfit_v2.1/v2.1')
addpath('minq5')


%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data
data.path = 'testdata/';  % path within EMdatapathname to find experimental patterns and to save data
imagenames = [
    "EDAX-Ni_Scan1_x126y0.png"
    "EDAX-Ni_Scan1_x126y2.png"
    "EDAX-Ni_Scan1_x126y3.png"
    ];  % name of image files to loop through (located in path)
data.binning = 8;        % binning mode (1, 2, 4, 8). Simulated image size will be numsx/binning x numsy/binning

% Define euler angles (each row phi1, PHI, phi2 in degrees)
euler = [
	35.2655 109.96 205.439
	35.2655 109.96 205.439
	35.2655 109.96 205.439
    ];

% Output data (outputs separate file for each image)
outputprefix = 'output';    % output file will be [outputprefix]_[imagename(ii) without extension]_v#.h5/txt

% Starting guesses for PC (usually from Hough data)
data.L = 15870.0234;    % in microns
data.xpc = 3.48576;      % in px
data.ypc = 114.20352;     % in px

% Define pseudosymmetry variants to check
data.pseudosym = [
    1 0 0 90;
    0 1 0 90
    ];     % [axis_x axis_y axis_z angle(deg)]


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
data.tmppath = 'tmp/';   % path within EMdatapathname for creation of temp files  (you need to manually create this directory)

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



N = length(imagenames);     % number of patterns to test
Np = size(data.pseudosym,1);     % number of pseudosymmetry variants


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
        
        % check pseudosymmetry variants
        for jj=1:Np
            h5output = strcat(outputprefix,'_',s,'_v',num2str(jj),'.h5');
            h5path = fullfile(data.homepath,data.path,h5output);
            if exist(h5path,'file')==2
                error('File already exists: %s',h5path);
            end
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
                
        % check pseudosymmetry variants
        for jj=1:Np
            logoutput = strcat(outputprefix,'_',s,'_v',num2str(jj),'.txt');
            logpath = fullfile(data.homepath,data.path,logoutput);
            if exist(logpath,'file')==2
                error('File already exists: %s',h5path);
            end
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
    data.outputprefix = strcat(outputprefix,'_',s);

    
    %%% Run PC fitting
    PCfit_snobfit_pseudo(data)
end


% read in data from each run
L = zeros(Np,N);     % initialize
xpc = zeros(Np,N);
ypc = zeros(Np,N);
dp = zeros(Np,N);

if data.h5outputlogical == 1
    for ii=1:N
        % figure out names of output files
        index = strfind(imagenames(ii),'.');
        s = char(imagenames(ii));
        s = s(1:index-1);
        h5output = strcat(outputprefix,'_',s,'.h5');
        
        % read in .h5 file
        h5path = fullfile(data.homepath,data.path,h5output);
        L(1,ii) = h5read(h5path,'/Data/L');
        xpc(1,ii) = h5read(h5path,'/Data/xpc');
        ypc(1,ii) = h5read(h5path,'/Data/ypc');
        dp(1,ii) = h5read(h5path,'/Data/dp');
        
        
        % pseudosymmetry variants
        for jj=1:Np
            h5output = strcat(outputprefix,'_',s,'_v',num2str(jj),'.h5');
            h5path = fullfile(data.homepath,data.path,h5output);
            L(jj+1,ii) = h5read(h5path,'/Data/L');
            xpc(jj+1,ii) = h5read(h5path,'/Data/xpc');
            ypc(jj+1,ii) = h5read(h5path,'/Data/ypc');
            dp(jj+1,ii) = h5read(h5path,'/Data/dp');
        end
    end
end


% figure out best pseudosym variant
[~,ibest] = max(dp);    % get indices of best variant for each pattern
L_best = zeros(1,N);    % initialize
xpc_best = zeros(1,N);
ypc_best = zeros(1,N);
for ii=1:N
    ix = ibest(ii);     % index of best variant for this pattern
    L_best(ii) = L(ix,ii);
    xpc_best(ii) = xpc(ix,ii);
    ypc_best(ii) = ypc(ix,ii);
end


% compute average PC
L_ave = sum(L_best)/N;
L_ave_frac = L_ave/(data.numsx*data.delta);
xpc_ave = sum(xpc_best)/N;
xpc_ave_frac = xpc_ave/data.numsx + 0.5;
ypc_ave = sum(ypc_best)/N;
ypc_ave_frac = ypc_ave/data.numsy + 0.5;


% print
fprintf('\n\nFitted pattern center (average of %g patterns):\n',N);
fprintf('  xpc = %.4f (X* = %.6f)\n',xpc_ave,xpc_ave_frac);
fprintf('  ypc = %.4f (Y* = %.6f)\n',ypc_ave,ypc_ave_frac);
fprintf('  L = %.2f (Z* = %.6f)\n\n',L_ave,L_ave_frac);


%%% write to file
fid = fopen(summarypath,'w');

fprintf(fid,'L:\n');
% original variant
fprintf(fid,'    v0:');
for ii=1:N  
    fprintf(fid,' %10.2f,',L(1,ii));
end
% pseudosym variants
for jj=1:Np
    fprintf(fid,'\n    v%g:',jj);
    for ii=1:N  
        fprintf(fid,' %10.2f,',L(jj+1,ii));
    end
end
% best variant
fprintf(fid,'\n  best:');
for ii=1:N  
    fprintf(fid,' %10.2f,',L_best(ii));
end

fprintf(fid,'\n\nxpc:\n');
% original variant
fprintf(fid,'    v0:');
for ii=1:N  
    fprintf(fid,' %10.4f,',xpc(1,ii));
end
% pseudosym variants
for jj=1:Np
    fprintf(fid,'\n    v%g:',jj);
    for ii=1:N  
        fprintf(fid,' %10.4f,',xpc(jj+1,ii));
    end
end
% best variant
fprintf(fid,'\n  best:');
for ii=1:N  
    fprintf(fid,' %10.4f,',xpc_best(ii));
end

fprintf(fid,'\n\nypc:\n');
% original variant
fprintf(fid,'    v0:');
for ii=1:N  
    fprintf(fid,' %10.4f,',ypc(1,ii));
end
% pseudosym variants
for jj=1:Np
    fprintf(fid,'\n    v%g:',jj);
    for ii=1:N  
        fprintf(fid,' %10.4f,',ypc(jj+1,ii));
    end
end
% best variant
fprintf(fid,'\n  best:');
for ii=1:N  
    fprintf(fid,' %10.4f,',ypc_best(ii));
end

fprintf(fid,'\n\ndp:\n');
% original variant
fprintf(fid,'    v0:');
for ii=1:N  
    fprintf(fid,' %10.6f,',dp(1,ii));
end
% pseudosym variants
for jj=1:Np
    fprintf(fid,'\n    v%g:',jj);
    for ii=1:N  
        fprintf(fid,' %10.6f,',dp(jj+1,ii));
    end
end

fprintf(fid,'\n\nFitted pattern center (average of %g patterns):\n',N);
fprintf(fid,'  xpc = %.4f (X* = %.6f)\n',xpc_ave,xpc_ave_frac);
fprintf(fid,'  ypc = %.4f (Y* = %.6f)\n',ypc_ave,ypc_ave_frac);
fprintf(fid,'  L = %.2f (Z* = %.6f)',L_ave,L_ave_frac);
fprintf(fid,'\n');
fclose(fid);

fprintf('Summarized data stored in txt file: %s\n',summarypath);


