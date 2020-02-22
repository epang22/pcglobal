%TESTERRORSTIMINGS
% Figure out parameters for FindParameters.m
% Test performance of EMEBSDDI, reporting errors and timings.
% Loop through:
% -nthreads: number of threads
% -numsingle: number of patterns arranged in column for dot product calcs on GPU
% -N_ori: number of orientations
% 2/21/20 (Edward Pang, MIT)


clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters to loop through (vectors)
nthreads = [6 8];        % number of threads
numsingle = [32 80];     % number of patterns arranged in column for dot product calcs on GPU (multiples of 16 perform better)
N_ori = [1 125];         % number of orientations

% run each this many times (certain combinations work only some of the time)
N = 2;

% Pattern size
data.numsx = 480;    % number of CCD pixels along x and y
data.numsy = 480;
data.binning = 8;

% Any master pattern on your computer
data.energymin = 10.0;   % energy range in the intensity summation (keV)
data.energymax = 20.0;
data.masterfile = 'DItutorial_Ni-master-20kV/Ni-master-20kV.h5';   % master pattern input file, path relative to EMdatapathname



%%% Parameters you don't need to change often %%%
% Paths for this computer
homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)
tmppath = 'tmp/';   % path within EMdatapathname for creation of temp files (you need to manually create this directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% clear tmp folder
[~,~] = system(sprintf('rm %s%s*',homepath,tmppath));     % delete all files


% Parameters you don't need to touch (values don't matter)
data.L = 15000;    % in microns
data.xpc = 0;         % in px
data.ypc = 80;       % in px
data.thetac = 5;   % tilt angle of the camera (positive below horizontal, degrees)
data.delta = 50;   % CCD pixel size on the scintillator surface (microns)
data.omega = 0;      % angle between normal of sample and detector
data.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
data.gammavalue = 0.33;  % gamma correction factor
data.maskpattern = 'y';        % use circular mask? y or n
data.hipassw = 0.05;    % hi pass filter w param; 0.05 is reasonable
data.nregions = 10;     % # of regions for adaptive histogram equalization
data.r = min(data.numsx,data.numsy)/(2*data.binning);        % radius of circular mask (after binning)
N_patterns = 1;		% Number of patterns to index


% create patterns and store in .data file
patterns = randi([0 255],data.numsy/data.binning,data.numsx/data.binning,N_patterns);
outname = 'test.data';
outpath = fullfile(homepath,tmppath,outname);
makebinary_multiimg(patterns,outpath);


% store some data for emebsddi_dp_multi3.m
data.homepath = homepath;
data.tmppath = tmppath;
data.img_exp = outname;
data.ipf_ht = 1;
data.ipf_wd = N_patterns;

% figure out some things
n_nthreads = length(nthreads);
n_numsingle = length(numsingle);
n_N_ori = length(N_ori);



% print header
fprintf('nthreads    N_ori  numsingle   # errors?  time(s)    ori/s\n');


% loop through parameters and call EMEBSDDI
for ii=1:n_nthreads
    data.nthreads = nthreads(ii);
    
    for jj=1:n_numsingle
        data.numdictsingle = numsingle(jj);   % number of dictionary files arranged in column for dp on GPU (multiples of 16 perform better)
        data.numexptsingle = numsingle(jj);   % number of experiment files "

        for kk=1:n_N_ori
            fprintf('%8.0f %8.0f %10.0f    ',nthreads(ii),numsingle(jj),N_ori(kk));
            
            % make N euler angles
            euler = zeros(N_ori(kk),3);
            euler(:,1) = linspace(0,45,N_ori(kk));
                
            % N trials of each parameter set
            timetotal = 0;  % keep track of total time
            errortally = 0;     % keep track of number of errors
            for ll=1:N
                try
                    t = tic;
                    dp = emebsddi_dp_multi3(euler,data,0);
                    time = toc(t);
                    timetotal = timetotal + time;
                    
                catch
                    errortally = errortally + 1;
                end
                
            end
            
            % compute average time
            if errortally<N
                timeave = timetotal/(N-errortally);     % average time of non-error runs
            else
                timeave = 0;    % if all runs give error
            end
            
            fprintf('%8.0f %8.2f %8.2f\n',errortally,timeave,N_ori(kk)/timeave);
        end
    end
end



% delete image .data file
[~,~] = system(sprintf('rm %s',outpath));




