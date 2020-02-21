function emebsd_singlepattern(L,xpc,ypc,euler,outpath,options)
%EMEBSD_WRAPPER_FUN
% Simulate a single pattern using EMsoft's EMEBSD program and save to disk
% 2/13/20
%
% Inputs: 
% -L (um)
% -xpc (px)
% -ypc (px)
% -euler=[phi1 PHI phi2] (deg)
% -outpath: full path of image name to output
% -options*
%
% *options is a struct containing the following fields:
% options.homepath = '/home/jonathan/epang/EMsoftfiles/EMdata/';  % path to EMdatapathname
% options.tmppath = 'tmp/';  % path within EMdatapathname for temp files and to find img_exp
% options.thetac = 10.0;   % tilt angle of the camera (positive below horizontal, degrees)
% options.delta = 59.2;   % CCD pixel size on the scintillator surface (microns)
% options.numsx = 480;    % number of CCD pixels along x and y
% options.numsy = 480;
% options.omega = 0;      % angle between normal of sample and detector
% options.energymin = 10.0;   % energy range in the intensity summation (keV)
% options.energymax = 20.0;
% options.eulerconvention = 'tsl';    % 'tsl' or 'hkl' Euler angle convention parameter
% options.masterfile = 'DItutorial_Ni-master-20kV/Ni-master-20kV.h5';   % master pattern input file, path relative to EMdatapathname
% options.poisson = 'n';      % include poisson noise? (y/n)
% options.binning = 8;        % binning mode (1, 2, 4, 8)
% options.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
% options.gammavalue = 0.33;  % gamma correction factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %


% some parameters that you don't need to specify
nthreads = 1;
maskpattern = 'n';  % apply circular mask to data? (y/n)
maskradius = 240;   % in pixels, AFTER binning operation


% error checking
if exist(outpath,'file')
    error('File already exists: %s',outpath);
end
if size(euler,1)>1
    error('More than 1 set of euler angles specified.');
end


% Name paths
eulerpath = fullfile(options.homepath,options.tmppath,'euler.txt');
eulerpathshort = fullfile(options.tmppath,'euler.txt');
nmlpath = fullfile(options.homepath,options.tmppath,'EMEBSD.nml');
h5path = fullfile(options.homepath,options.tmppath,'EBSDout.h5');
h5pathshort = fullfile(options.tmppath,'EBSDout.h5');


% Create euler.txt file
fid = fopen(eulerpath,'w');
fprintf(fid,'eu\n');
fprintf(fid,'1\n');
fprintf(fid,' %g %g %g\n',euler(1),euler(2),euler(3));
fclose(fid);


% Create EMEBSD.nml file
fid = fopen(nmlpath,'w');
fprintf(fid,'&EBSDdata\n');
fprintf(fid,' L = %.2f,\n',L);
fprintf(fid,' thetac = %.1f,\n',options.thetac);
fprintf(fid,' delta = %.1f,\n',options.delta);
fprintf(fid,' numsx = %g,\n',options.numsx);
fprintf(fid,' numsy = %g,\n',options.numsy);
fprintf(fid,' xpc = %.4f,\n',xpc);
fprintf(fid,' ypc = %.4f,\n',ypc);
fprintf(fid,' omega = %.1f,\n',options.omega);
fprintf(fid,' alphaBD = 0.0,\n');
fprintf(fid,' energymin = %.1f,\n',options.energymin);
fprintf(fid,' energymax = %.1f,\n',options.energymax);
fprintf(fid,' includebackground = ''y'',\n');
fprintf(fid,' anglefile = ''%s'',\n',eulerpathshort);
fprintf(fid,' anglefiletype = ''orientations'',\n');
fprintf(fid,' eulerconvention = ''%s'',\n',options.eulerconvention);
fprintf(fid,' masterfile = ''%s'',\n',options.masterfile);
fprintf(fid,' energyfile = ''%s'',\n',options.masterfile);
fprintf(fid,' datafile = ''%s'',\n',h5pathshort);
fprintf(fid,' bitdepth = ''8bit'',\n');
fprintf(fid,' beamcurrent = 150.0,\n');
fprintf(fid,' dwelltime = 100.0,\n');
fprintf(fid,' poisson = ''%s'',\n',options.poisson);
fprintf(fid,' binning = %g,\n',options.binning);
fprintf(fid,' applyDeformation = ''n'',\n');
fprintf(fid,' Ftensor = 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0,\n');
fprintf(fid,' scalingmode = ''%s'',\n',options.scalingmode);
fprintf(fid,' gammavalue = %g,\n',options.gammavalue);
fprintf(fid,' makedictionary = ''n'',\n');
fprintf(fid,' maskpattern = ''%s'',\n',maskpattern);
fprintf(fid,' maskradius = %g,\n',maskradius);
fprintf(fid,' hipassw = 0.05,\n');
fprintf(fid,' nregions = 10,\n');
fprintf(fid,' nthreads = %g,\n',nthreads);
fprintf(fid,' /\n');
fclose(fid);


% Run EMsoft EMEBSD
[~,~] = system(sprintf('EMEBSD %s',nmlpath));


% Extract images and save to file
data = h5read(h5path,'/EMData/EBSD/EBSDPatterns');

% % use this line if EMsoft v4.0
% img = flipud(data');     % invert and flip array so that image is oriented properly (TSL convention)
% use this line if EMsoft v4.3
img = data';

imwrite(img,outpath);



% delete files
[~,~] = system(sprintf('rm %s',eulerpath));
[~,~] = system(sprintf('rm %s',nmlpath));
[~,~] = system(sprintf('rm %s',h5path));

