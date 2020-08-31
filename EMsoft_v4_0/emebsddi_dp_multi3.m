function dp = emebsddi_dp_multi3(euler,options,verbose)
%EMEBSDDI_DP_MULTI Compute dot product by calling EMsoft EMEBSDDI program
% Run EMsoft EMEBSDDI program to simulate EBSD pattern for multiple orientations,
% single set of detector parameters, compute dot product with a multiple
% experimental patterns
% Compared to _multi2, can accommodate maps with arbitrary ipf
% Inputs: 
% >Np=number of patterns stored in options.img_exp2
% >euler=[phi1 PHI phi2
%         phi1 PHI phi2
%         ...          ] (deg)
% >options*
% >verbose: 1=print emsoft output, 0=no (optional, will print output if not specified)
%
% *options is a struct containing the following fields:
% options.homepath = '/home/jonathan/epang/EMsoftfiles/EMdata/';  % path to EMdatapathname
% options.tmppath = 'tmp/';  % path within EMdatapathname for temp files
% options.img_exp = 'test.data';    % name of .data file containing experimental patterns (located in tmppath)
% options.ipf_ht = 150;     % height of data set in pattern input file
% options.ipf_wd = 150;     % width
% options.L = 15594.81;    % in microns
% options.xpc = 0.9558;         % in px
% options.ypc = 65.8671;       % in px
% options.thetac = 10.0;   % title angle of the camera (positive below horizontal, degrees)
% options.delta = 59.2;   % CCD pixel size on the scintillator surface (microns)
% options.numsx = 480;    % number of CCD pixels along x and y
% options.numsy = 480;
% options.omega = 0;      % angle between normal of sample and detector
% options.energymin = 10.0;   % energy range in the intensity summation (keV)
% options.energymax = 20.0;
% options.masterfile = 'DItutorial_Ni-master-20kV/Ni-master-20kV.h5';   % master pattern input file, path relative to EMdatapathname
% options.binning = 8;        % binning mode (1, 2, 4, 8)
% options.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
% options.gammavalue = 0.33;  % gamma correction factor
% options.nthreads = 8;
% options.maskpattern = 'y';        % use circular mask? y or n
% options.r = 30;        % radius of circular mask for computing dp
% options.hipassw = 0.05;        % hi pass filter w param; 0.05 is reasonable
% options.nregions = 10;      % # of regions for adaptive histogram equalization
% options.numdictsingle = 1024;   % number of dictionary files arranged in column for dp on GPU (multiples of 16 perform better)
% options.numexptsingle = 1024;   % number of experiment files "
%
% Outputs:
% >dp: N_angles x Np array of dp in same order as euler, each col is a map point
%
% Original: 1/12/20 (Edward Pang, MIT)
% Change log:
% 8/31/20 ELP: fprintf change %g (6 sig figs) to %.8f for euler angles, L, xpc, ypc


% Name paths
eulerpath = fullfile(options.homepath,options.tmppath,'euler.txt');
eulerpathshort = fullfile(options.tmppath,'euler.txt');
nmlpath = fullfile(options.homepath,options.tmppath,'EMEBSDDI.nml');
h5path = fullfile(options.homepath,options.tmppath,'DIoutput.h5');
h5pathshort = fullfile(options.tmppath,'DIoutput.h5');
outpathshort = fullfile(options.tmppath,options.img_exp);


% Figure out how many angles
N_angles = size(euler,1);   


% Create EMEBSDDI.nml file
fid = fopen(nmlpath,'w');
fprintf(fid,' &EBSDIndexingdata\n');
fprintf(fid,' indexingmode = ''dynamic'',\n');
fprintf(fid,' Notify = ''off'',\n');
fprintf(fid,' ipf_ht = %g,\n',options.ipf_ht);    % height of data set in pattern input file
fprintf(fid,' ipf_wd = %g,\n',options.ipf_wd);          % width
fprintf(fid,' ROI = 0 0 0 0,\n');       % leave all at 0 for full field of view
fprintf(fid,' stepX = 1.0,\n');         % X and Y sampling step sizes
fprintf(fid,' stepY = 1.0,\n');
fprintf(fid,' nnk = %g,\n',N_angles);          % # of top dot products to save (this defines the max # of orientations)
fprintf(fid,' nnav = 1,\n');
fprintf(fid,' nosm = 1,\n');
fprintf(fid,' maskfile = ''undefined'',\n');
fprintf(fid,' maskpattern = ''%s'',\n',options.maskpattern);
fprintf(fid,' maskradius = %.0f,\n',options.r);
fprintf(fid,' hipassw = %g,\n',options.hipassw);
fprintf(fid,' nregions = %g,\n',options.nregions);

fprintf(fid,' ncubochoric = 40,\n');
fprintf(fid,' L = %.8f,\n',options.L);
fprintf(fid,' thetac = %g,\n',options.thetac);
fprintf(fid,' delta = %g,\n',options.delta);
fprintf(fid,' numsx = %g,\n',options.numsx);
fprintf(fid,' numsy = %g,\n',options.numsy);
fprintf(fid,' xpc = %.8f,\n',options.xpc);
fprintf(fid,' ypc = %.8f,\n',options.ypc);
fprintf(fid,' omega = %g,\n',options.omega);
fprintf(fid,' energymin = %g,\n',options.energymin);
fprintf(fid,' energymax = %g,\n',options.energymax);
fprintf(fid,' energyaverage = 0,\n');
fprintf(fid,' spatialaverage = ''n'',\n');
fprintf(fid,' beamcurrent = 150.0,\n');
fprintf(fid,' dwelltime = 100.0,\n');
fprintf(fid,' binning = %g,\n',options.binning);
fprintf(fid,' scalingmode = ''%s'',\n',options.scalingmode);
fprintf(fid,' gammavalue = %g,\n',options.gammavalue);

fprintf(fid,' exptfile = ''%s'',\n',outpathshort);
fprintf(fid,' inputtype = ''Binary'',\n');
fprintf(fid,' HDFstrings = '''' '''' '''' '''' '''' '''' '''' '''' '''' '''',\n');

fprintf(fid,' tmpfile = ''EMEBSDDict_tmp.data'',\n');
fprintf(fid,' keeptmpfile = ''n'',\n');
fprintf(fid,' datafile = ''%s'',\n',h5pathshort);
fprintf(fid,' ctffile = ''undefined'',\n');
fprintf(fid,' avctffile = ''undefined'',\n');
fprintf(fid,' eulerfile = ''%s'',\n',eulerpathshort);

fprintf(fid,' dictfile = ''undefined'',\n');

fprintf(fid,' masterfile = ''%s'',\n',options.masterfile);

fprintf(fid,' numdictsingle = %g,\n',options.numdictsingle);
fprintf(fid,' numexptsingle = %g,\n',options.numexptsingle);
fprintf(fid,' nthreads = %g,\n',options.nthreads);
fprintf(fid,' platid = 1,\n');
fprintf(fid,' devid = 1,\n');
fprintf(fid,' multidevid = 1 0 0 0 0 0 0 0,\n');
fprintf(fid,' usenumd = 1,\n');

fprintf(fid,' /\n');
fclose(fid);


% Create euler.txt file
fid = fopen(eulerpath,'w');
fprintf(fid,'eu\n');
fprintf(fid,'%.0f\n',size(euler,1));
for ii=1:N_angles
    fprintf(fid,' %.8f %.8f %.8f\n',euler(ii,1),euler(ii,2),euler(ii,3));
end
fclose(fid);


% Run EMsoft EMEBSDDI
if nargin<3
    system(sprintf('EMEBSDDI %s',nmlpath));     % verbose
else
    if verbose==1
        system(sprintf('EMEBSDDI %s',nmlpath));     % verbose
    else
        [~,~] = system(sprintf('EMEBSDDI %s',nmlpath));     % don't print emsoft output
    end
end


% Read in data
% sometimes can't find h5path immediately after running EMEBSDDI, these lines of code try to remedy that
t = tic;    
while exist(h5path,'file')==0 && toc(t)<1
    pause(0.1);     % pause 0.1 seconds and try again
end
if exist(h5path,'file')==2
    TopMatchIndices = h5read(h5path,'/Scan 1/EBSD/Data/TopMatchIndices');
    TopDotProductList = h5read(h5path,'/Scan 1/EBSD/Data/TopDotProductList');
    NumExptPatterns = h5read(h5path,'/Scan 1/EBSD/Data/NumExptPatterns');
else
    error('h5path not found!');
end

% Chop empty data
if size(TopMatchIndices,2) > NumExptPatterns
    TopMatchIndices(:,NumExptPatterns+1:end) = [];
    TopDotProductList(:,NumExptPatterns+1:end) = [];
end

% Resort dot products in same order as EulerAngles
dp = zeros(N_angles,NumExptPatterns);  % initialize
for ii=1:NumExptPatterns
    for jj=1:N_angles
        index = TopMatchIndices(jj,ii);
        if index>N_angles
            error('index > N_angles. index = %g',index)
        end
        dp(index,ii) = TopDotProductList(jj,ii);
    end
end
   


% delete input files
[~,~] = system(sprintf('rm %s',eulerpath));
[~,~] = system(sprintf('rm %s',nmlpath));
[~,~] = system(sprintf('rm %s',h5path));

