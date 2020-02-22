%FINDPARAMETERS
% Linescans through PC/orientation space to determine scaling factors and 
% min step sizes for RunPCfit
% Generates a simulated pattern with specified parameters, rather than 
% loading in an experimental pattern, to ensure that you scan through the 
% peak tip to get accurate FWHM
% 2/13/20

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Center of linescan
%  Input any values for pattern center in the vicinity you care about
%  Euler angles don't really matter
L = 15000;  % in microns
xpc = 0;    % in px
ypc = 80;   % inpx
euler = [283.164 93.5072 339.491];  % euler angles (deg)

% linescan parameters
N = 50;     % # of points in each linescan
delta_L = 2500;  % linescan will check +/- this amount, in microns
delta_xpc = 25;   % " in px
delta_ypc = 25;   % " in px
delta_angles = 5;     % " in degrees

% EMsoft parameters for the master pattern
energymin = 10.0;   % energy range in the intensity summation (keV)
energymax = 20.0;
masterfile = 'DItutorial_Ni-master-20kV/Ni-master-20kV.h5';   % master pattern input file, path relative to EMdatapathname

%%% EMEBSDDI parameters
% These parameters affect how long the calculations will take. Also, some combinations give errors (bugs in EMsoft).
% Use TestErrorsTimings.m to determine appropriate parameters
nthreads = 8; % number of threads
% number of files arranged in column for dp on GPU (multiples of 16 perform better)
numsingle1 = 32;   % Must be compatible with nthreads, N_ori=1
numsingle2 = 80;   % Must be compatible with nthreads, N_ori=3*N

progress = 10;  % print progress every this many iterations (out of 3*N)



%%% Parameters you don't need to change often %%%
% Paths for this computer
data.homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)
data.tmppath = 'tmp/';   % path within EMdatapathname for creation of temp files (you need to manually create this directory)

% Detector parameters
data.thetac = 10;   % tilt angle of the camera (positive below horizontal, degrees)
data.delta = 59.2;   % CCD pixel size on the scintillator surface (microns)
data.numsx = 480;    % number of CCD pixels along x and y
data.numsy = 480;
data.omega = 0;      % angle between normal of sample and detector
binning = 1;        % binning mode (1, 2, 4, 8). Simulated image size will be numsx/binning x numsy/binning

% Some EMsoft parameters for dot product computations
eulerconvention = 'tsl';    % 'tsl' or 'hkl' Euler angle convention parameter
poisson = 'n';      % include poisson noise? (y/n)
scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
gammavalue = 0.33;  % gamma correction factor
maskpattern = 'y';        % use circular mask? y or n
r = min(numsx,numsy)/(2*binning);        % radius of circular mask (after binning)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %



totaltime = tic;



% clear tmp folder
[~,~] = system(sprintf('rm %s%s*',homepath,tmppath));     % delete all files


% store info in struct
options.homepath = homepath;
options.tmppath = tmppath;
outname = 'test.data';
options.img_exp = outname;
options.thetac = thetac;
options.delta = delta;
options.numsx = numsx;
options.numsy = numsy;
options.omega = omega;
options.energymin = energymin;
options.energymax = energymax;
options.masterfile = masterfile;
options.eulerconvention = eulerconvention;
options.poisson = poisson;
options.binning = binning;
options.scalingmode = scalingmode;
options.gammavalue = gammavalue;
options.nthreads = nthreads;
options.maskpattern = maskpattern;
options.r = r;


% Simulate pattern using EMEBSD
imgpath = fullfile(homepath,tmppath,'test.png');
emebsd_singlepattern(L,xpc,ypc,euler,imgpath,options);


% Convert simulated pattern to binary .data file
outpath = fullfile(homepath,tmppath,outname);
makebinary_singleimg(imgpath,outpath,options)
    

fprintf('Initialization complete.\n');



%%% compute Rx,Ry,Rz linescans
fprintf('Computing orientation linescans...\n');

% update some parameters for EMEBSDDI
options.numdictsingle = numsingle2;
options.numexptsingle = numsingle2;
options.nnk = 3*N;

% R values for linescan
test_R = linspace(-tand(delta_angles/2),tand(delta_angles/2),N);  

% compute linescan euler angles
test_angles_Rx = createorigrid_1d(test_R,1,euler);   % Rx linescan angles
test_angles_Ry = createorigrid_1d(test_R,2,euler);   % Rx linescan angles
test_angles_Rz = createorigrid_1d(test_R,3,euler);   % Rx linescan angles
test_angles_all = [test_angles_Rx; test_angles_Ry; test_angles_Rz];     % concatenate, to feed all into EMEBSDDI at same time

% compute dot products
dp_R_all = emebsddi_dp(L,xpc,ypc,test_angles_all,options);    

% split into Rx, Ry, Rz
dp_Rx = dp_R_all(1:N);
dp_Ry = dp_R_all(N+1:2*N);
dp_Rz = dp_R_all(2*N+1:3*N);



%%% Compute PC linescans
fprintf('Computing PC linescans...\n');

% update some parameters for EMEBSDDI
options.numdictsingle = numsingle1;
options.numexptsingle = numsingle1;
options.nnk = 1;

looptime = tic;     % start timer for loop
% Scan in L
test_L = linspace(L-delta_L,L+delta_L,N);    % L values of linescan
dp_L = zeros(1,N);     % initialize array to store dp
for ii=1:N
    dp_L(ii) = emebsddi_dp(test_L(ii),xpc,ypc,euler,options);  % compute dp
    displayprogress(ii,3*N,toc(looptime),progress);   % print progress
end
    
% scan in xpc
test_xpc = linspace(xpc-delta_xpc,xpc+delta_xpc,N);    % xpc values of linescan
dp_xpc = zeros(1,N);     % initialize array to store dp
for ii=1:N
    dp_xpc(ii) = emebsddi_dp(L,test_xpc(ii),ypc,euler,options);
    displayprogress(N+ii,3*N,toc(looptime),progress);   % print progress
end
    
% scan in ypc
test_ypc = linspace(ypc-delta_ypc,ypc+delta_ypc,N);    % ypc values of linescan
dp_ypc = zeros(1,N);     % initialize array to store dp
for ii=1:N
    dp_ypc(ii) = emebsddi_dp(L,xpc,test_ypc(ii),euler,options);
    displayprogress(2*N+ii,3*N,toc(looptime),progress);   % print progress
end
    


%%% Compute FWHM info
[width_xpc, halfmax_xpc, center_xpc] = fwhm(test_xpc,dp_xpc);
[width_ypc, halfmax_ypc, center_ypc] = fwhm(test_ypc,dp_ypc);
[width_L, halfmax_L, center_L] = fwhm(test_L,dp_L);
[width_Rx, halfmax_Rx, center_Rx] = fwhm(test_R,dp_Rx);
[width_Ry, halfmax_Ry, center_Ry] = fwhm(test_R,dp_Ry);
[width_Rz, halfmax_Rz, center_Rz] = fwhm(test_R,dp_Rz);
width_R = (width_Rx+width_Ry+width_Rz)/3;   % compute average
width_R_deg = 4*atand(width_R/2);
halfmax_R = (halfmax_Rx+halfmax_Ry+halfmax_Rz)/3;   % compute average
center_R = (center_Rx+center_Ry+center_Rz)/3;   % compute average


%%% Compute recommended parameters
% just scale from our paper
width_L2 = (width_L/1200)*225;
width_xpc2 = (width_xpc/12)*2.25;
width_ypc2 = (width_ypc/12)*2.25;
width_R_deg2 = (width_R_deg/2.4)*0.43;
Lstep_min = (width_L/1200)*0.25;
xpcstep_min = (width_xpc/12)*0.0025;
ypcstep_min = (width_ypc/12)*0.0025;
anglestep_min = (width_R_deg/2.4)*0.0005;


% print
fprintf('\n');
fprintf('Use the following parameters for RunPCfit:\n');
fprintf('  Lstep = %.2f\n',width_L2);
fprintf('  xpcstep = %.2f\n',width_xpc2);
fprintf('  ypcstep = %.2f\n',width_ypc2);
fprintf('  anglestep = %.2f\n',width_R_deg2);
fprintf('  Lstep_min = %.3f\n',Lstep_min);
fprintf('  xpcstep_min = %.4f\n',xpcstep_min);
fprintf('  ypcstep_min = %.4f\n',ypcstep_min);
fprintf('  anglestep_min = %.5f\n',anglestep_min);
fprintf('\n');

    
%%% Plots
% plot parameters
linewidth = 1.5;
axislinewidth = 1.2;
fontaxislabel = 14;
fontticklabel = 11;

% L
figure;
plot(test_L,dp_L,'LineWidth',linewidth);
hold on
plot([center_L-width_L/2 center_L+width_L/2],halfmax_L*[1 1],'-k','LineWidth',linewidth);
set(gca,'FontSize',fontticklabel,'LineWidth',axislinewidth);
xlabel('L (μm)','FontSize',fontaxislabel);
ylabel('NDP','FontSize',fontaxislabel);

% xpc
figure;
plot(test_xpc,dp_xpc,'LineWidth',linewidth);
hold on
plot([center_xpc-width_xpc/2 center_xpc+width_xpc/2],halfmax_xpc*[1 1],'-k','LineWidth',linewidth);
set(gca,'FontSize',fontticklabel,'LineWidth',axislinewidth);
xlabel('xpc (px)','FontSize',fontaxislabel);
ylabel('NDP','FontSize',fontaxislabel);

% ypc
figure;
plot(test_ypc,dp_ypc,'LineWidth',linewidth);
hold on
plot([center_ypc-width_ypc/2 center_ypc+width_ypc/2],halfmax_ypc*[1 1],'-k','LineWidth',linewidth);
set(gca,'FontSize',fontticklabel,'LineWidth',axislinewidth);
xlabel('ypc (px)','FontSize',fontaxislabel);
ylabel('NDP','FontSize',fontaxislabel);

% orientation
figure;
plot(test_R,dp_Rx,'LineWidth',linewidth);
hold on
plot(test_R,dp_Ry,'LineWidth',linewidth);
plot(test_R,dp_Rz,'LineWidth',linewidth);
plot([center_R-width_R/2 center_R+width_R/2],halfmax_R*[1 1],'-k','LineWidth',linewidth);
set(gca,'FontSize',fontticklabel,'LineWidth',axislinewidth);
xlabel('R','FontSize',fontaxislabel);
ylabel('NDP','FontSize',fontaxislabel);

    % add misorientation second axis
    set(gca,'OuterPosition',[0 0 1 0.95]);  % make sure top axis doesn't get cut off
    x1 = xlim;  % get current axis limits for calculating misorientation axis
    ax1 = gca; ax1_pos = ax1.Position;
    set(ax1,'box','off');    % remove upper x ticks on ax1
    ax1b = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
        'xtick',[],'ytick',[],'color','none');
    ax2 = axes('Position',ax1_pos,'Color','none','XAxisLocation','top','YAxisLocation','right',...
        'FontSize',fontticklabel,'LineWidth',axislinewidth);
    ax2.YAxis.Visible = 'off';    % define top axis
    xlim(ax2,[2*atand(x1(1)) 2*atand(x1(2))]);  % set limits misorientation in deg to match Rx
    xlabel(ax2,'Misorientation (deg)');




% delete image .data file
[~,~] = system(sprintf('rm %s',outpath));




toc(totaltime)




%%% Define function to display progress
function displayprogress(ii,N,time,progress)
% print progress, estimated time remaining to screen
% ii=iterations completed
% N=total iterations
% time=time in sec to complete ii iterations
% progress=only display progress if ii is a multiple of this number

if mod(ii,progress)==0
    percentcomplete = ii/N;
    esttimeremain = time*(N-ii)/ii;
    if esttimeremain<60
        fprintf('    %.1f%% completed. Estimated time remaining: %.1f seconds...\n',100*percentcomplete,esttimeremain);
    elseif esttimeremain<3600
        fprintf('    %.1f%% completed. Estimated time remaining: %.1f minutes...\n',100*percentcomplete,esttimeremain/60);
    else
        fprintf('    %.1f%% completed. Estimated time remaining: %.1f hours...\n',100*percentcomplete,esttimeremain/3600);
    end
end
        
end



%%% Define function to compute orientaiton
function test_angles = createorigrid_1d(R,dir,euler)
% Compute 1D grid of orientations
% Inputs:
% -R: vector of Rodrigues components
% -dir: grid direction? Rx=1, Ry=2, Rz=3
% -euler: orientation of grid center, in deg
% Output:
% -test_angles: Nx3 array of euler angles (deg)

N = length(R);
test_angles = zeros(N,3);    % Initialize array, contains output euler angles
q1 = eu2qu(euler*pi/180);  % convert euler angles of grid center to quaternions

% loop through each rodrigues component, convert to euler angles
for ss=1:N
    x = 0; y = 0; z = 0;
    if dir==1
        x = R(ss);    % Rodrigues vector components (grid centered about origin)
    elseif dir==2
        y = R(ss);    % Rodrigues vector components (grid centered about origin)
    elseif dir==3
        z = R(ss);    % Rodrigues vector components (grid centered about origin)
    else
        error('Unrecognized direction. dir must be 1, 2, or 3.');
    end
    
    w = 1/sqrt(1+x^2+y^2+z^2);          % Calculate first quaternion entry
    q2 = [w w*x w*y w*z];   % Quaternion of grid point (centered around north pole)

    qrot = [q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)-q1(4)*q2(4), ...
        q1(1)*q2(2)+q1(2)*q2(1)+q1(3)*q2(4)-q1(4)*q2(3), ...
        q1(1)*q2(3)+q1(3)*q2(1)-q1(2)*q2(4)+q1(4)*q2(2), ...
        q1(1)*q2(4)+q1(4)*q2(1)+q1(2)*q2(3)-q1(3)*q2(2)];   % quaternion grid rotated to be centered on starting Euler angles
    euler1 = qu2eu(qrot)*180/pi;     % Euler angles of grid point (in degrees)
    test_angles(ss,:) = euler1;
end 

end



%%% Define function to compute FWHM
function [width, halfmax, center] = fwhm(x,y)
% compute FWHM of a function (peak pointing up)
% Inputs: vectors x,y containing signal
% Outputs: 
% -width: width at FWHM
% -halfmax: y-value of FWHM
% -center: x-value at center of FWHM

[y_max, imax] = max(y);
y_min = min(y);
halfmax = (y_max+y_min)/2;

% find left side of FWHM
i = imax;
while sign(y(i)-halfmax)==sign(y(i+1)-halfmax)
    i = i-1;    % start at max point, move left until reach halfmax
end
ileft = i;  % index of left side of FWHM

% find right side of FWHM
i = imax;
while sign(y(i)-halfmax)==sign(y(i+1)-halfmax)
    i = i+1;    % start at max point, move left until reach halfmax
end
iright = i;

% compute FWHM
width = x(iright)-x(ileft);

% compute center of FWHM
center = (x(iright)+x(ileft))/2;

end

