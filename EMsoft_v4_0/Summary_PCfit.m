%SUMMARY_PCFIT
% Use this script to compute average PC from multiple PCfit runs. Also
% outputs a summary .txt file and diagnostic plots.
% Fill in INPUT PARAMETERS section with desired parameters
% 5/8/21 (Edward Pang, MIT)

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data
path = 'testdata/';  % path within EMdatapathname to find experimental patterns and to save data

% for a non-pseudosymmetric material:
h5names = [
    "output_EDAX-Ni_Scan1_x126y0.h5"
    "output_EDAX-Ni_Scan1_x126y2.h5"
	"output_EDAX-Ni_Scan1_x126y3.h5"
    ];  % list of .h5 files outputted by PCfit (located in path)

% % for a pseudosymmetric material:
% h5names = [
%     "output_EDAX-Ni_Scan1_x126y0.h5", "output_EDAX-Ni_Scan1_x126y0_v1.h5", "output_EDAX-Ni_Scan1_x126y0_v2.h5"
%     "output_EDAX-Ni_Scan1_x126y2.h5", "output_EDAX-Ni_Scan1_x126y2_v1.h5", "output_EDAX-Ni_Scan1_x126y2_v2.h5"
%     "output_EDAX-Ni_Scan1_x126y3.h5", "output_EDAX-Ni_Scan1_x126y3_v1.h5", "output_EDAX-Ni_Scan1_x126y3_v2.h5"
%     ];  % list of .h5 files outputted by PCfit (located in path)
%         % each row corresponds to a pattern, columns are pseudosym variants
        
        
% name of summary .txt output file (will be saved in path)
outname = 'output_summary_test.txt';
txtoutputlogical = 1; % output summary txt file? 1=yes, 0=no

%%% Parameters you don't need to change often %%%
% Paths for this computer
homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)

% Detector parameters (for computing X*,Y*,Z* in EDAX convention)
thetac = 10;   % tilt angle of the camera (positive below horizontal, degrees)
delta = 59.2;   % CCD pixel size on the scintillator surface (microns)
numsx = 480;    % number of CCD pixels along x and y
numsy = 480;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %


% name output path, check if file already exists
if txtoutputlogical==1
    summarypath = fullfile(homepath,path,outname);
    if exist(summarypath,'file')==2
        error('File already exists: %s',summarypath);
    end
end


% figure out number of files to read in
N = size(h5names,1);     % number of patterns to test
Np = size(h5names,2);     % number of pseudosym variants


% initialize some arrays
L = zeros(Np,N);     % PC info for all pseudosym variants
xpc = zeros(Np,N);
ypc = zeros(Np,N);
dp = zeros(Np,N);
L_best = zeros(1,N);    % PC info for best variant
xpc_best = zeros(1,N);
ypc_best = zeros(1,N);


% read in data from PCfit
for ii=1:N
    for jj=1:Np
        h5path = fullfile(homepath,path,h5names(ii,jj));
        L(jj,ii) = h5read(h5path,'/Data/L');
        xpc(jj,ii) = h5read(h5path,'/Data/xpc');
        ypc(jj,ii) = h5read(h5path,'/Data/ypc');
        dp(jj,ii) = h5read(h5path,'/Data/dp');
    end
end


% extract PC for best variant (highest dp)
if Np>1     % pseudosymmetric material
    [~,ibest] = max(dp);    % get indices of best variant for each pattern
    for ii=1:N
        ix = ibest(ii);     % index of best variant for this pattern
        L_best(ii) = L(ix,ii);
        xpc_best(ii) = xpc(ix,ii);
        ypc_best(ii) = ypc(ix,ii);
    end
else        % non-pseudosymmetric material
    % only one variant, just copy data
    L_best = L;
    xpc_best = xpc;
    ypc_best = ypc;
end


% compute average PC
L_ave = sum(L_best)/N;
L_ave_frac = L_ave/(numsx*delta);
xpc_ave = sum(xpc_best)/N;
xpc_ave_frac = xpc_ave/numsx + 0.5;
ypc_ave = sum(ypc_best)/N;
ypc_ave_frac = ypc_ave/numsy + 0.5;


% print
fprintf('\n\nFitted pattern center (average of %g patterns):\n',N);
fprintf('  xpc = %.4f (X* = %.6f)\n',xpc_ave,xpc_ave_frac);
fprintf('  ypc = %.4f (Y* = %.6f)\n',ypc_ave,ypc_ave_frac);
fprintf('  L = %.2f (Z* = %.6f)\n\n',L_ave,L_ave_frac);


% make plots
figure('pos',[300 200 1200 400]);

subplot(1,3,1);
plot(zeros(size(xpc_best)),xpc_best,'o');  % individual points
hold on
plot([-1 1],[xpc_ave xpc_ave],'--k');
ylabel('xpc (px)');

set(gca,'FontSize',12,'ticklength',[0.03 0.05],'XTick',[],'YMinorTick','on');

    % add fraction detector width second axis
    set(gca,'OuterPosition',[0 0 0.22 1]);  % make sure right axis doesn't get cut off
    yl = ylim;  % get current axis limits for calculating second axis
    ax1 = gca; ax1_pos = ax1.Position;
    set(ax1,'box','off');    % remove right y ticks on ax1
    ax2 = axes('Position',ax1_pos,'Color','none','XAxisLocation','top',...
        'YAxisLocation','right','FontSize',12,'XTick',[],'YMinorTick','on',...
        'ticklength',[0.03 0.05]);    % define top/right axis
    ylim(ax2,yl/numsx+0.5);  % set limits to match left axis
    ylabel(ax2,'X^* (EDAX)','FontSize',12,'Color','k');


subplot(1,3,2);
plot(zeros(size(ypc_best)),ypc_best,'o');  % individual points
hold on
plot([-1 1],[ypc_ave ypc_ave],'--k');
ylabel('ypc (px)');

set(gca,'FontSize',12,'ticklength',[0.03 0.05],'XTick',[],'YMinorTick','on');

    % add fraction detector width second axis
    set(gca,'OuterPosition',[0.33 0 0.22 1]);  % make sure right axis doesn't get cut off
    yl = ylim;  % get current axis limits for calculating second axis
    ax1 = gca; ax1_pos = ax1.Position;
    set(ax1,'box','off');    % remove right y ticks on ax1
    ax2 = axes('Position',ax1_pos,'Color','none','XAxisLocation','top',...
        'YAxisLocation','right','FontSize',12,'XTick',[],'YMinorTick','on',...
        'ticklength',[0.03 0.05]);    % define top/right axis
    ylim(ax2,yl/numsy+0.5);  % set limits to match left axis
    ylabel(ax2,'Y^* (EDAX)','FontSize',12,'Color','k');


subplot(1,3,3);
plot(zeros(size(L_best)),L_best,'o');  % individual points
hold on
plot([-1 1],[L_ave L_ave],'--k');
ylabel('L (μm)');

set(gca,'FontSize',12,'ticklength',[0.03 0.05],'XTick',[],'YMinorTick','on');
ax = gca; ax.YRuler.Exponent = 0;

    % add fraction detector width second axis
    set(gca,'OuterPosition',[0.66 0 0.26 1]);  % make sure right axis doesn't get cut off
    yl = ylim;  % get current axis limits for calculating second axis
    ax1 = gca; ax1_pos = ax1.Position;
    set(ax1,'box','off');    % remove right y ticks on ax1
    ax2 = axes('Position',ax1_pos,'Color','none','XAxisLocation','top',...
        'YAxisLocation','right','FontSize',12,'XTick',[],'YMinorTick','on',...
        'ticklength',[0.03 0.05]);    % define top/right axis
    ylim(ax2,yl/(numsx*delta));  % set limits to match left axis
    ylabel(ax2,'Z^* (EDAX)','FontSize',12,'Color','k');



%%% write to file (if specified)
if txtoutputlogical==1
    fid = fopen(summarypath,'w');

    if Np>1     % pseudosymmetric material
        fprintf(fid,'L:');
        for jj=1:Np
            fprintf(fid,'\n    v%g:',jj-1);
            for ii=1:N  
                fprintf(fid,' %10.2f,',L(jj,ii));
            end
        end
        % best variant
        fprintf(fid,'\n  best:');
        for ii=1:N  
            fprintf(fid,' %10.2f,',L_best(ii));
        end

        fprintf(fid,'\n\nxpc:');
        for jj=1:Np
            fprintf(fid,'\n    v%g:',jj-1);
            for ii=1:N  
                fprintf(fid,' %10.4f,',xpc(jj,ii));
            end
        end
        % best variant
        fprintf(fid,'\n  best:');
        for ii=1:N  
            fprintf(fid,' %10.4f,',xpc_best(ii));
        end

        fprintf(fid,'\n\nypc:');
        for jj=1:Np
            fprintf(fid,'\n    v%g:',jj-1);
            for ii=1:N  
                fprintf(fid,' %10.4f,',ypc(jj,ii));
            end
        end
        % best variant
        fprintf(fid,'\n  best:');
        for ii=1:N  
            fprintf(fid,' %10.4f,',ypc_best(ii));
        end

        fprintf(fid,'\n\ndp:');
        for jj=1:Np
            fprintf(fid,'\n    v%g:',jj-1);
            for ii=1:N  
                fprintf(fid,' %10.6f,',dp(jj,ii));
            end
        end

    else    % non-pseudosymmetric material
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

    end


    fprintf(fid,'\n\nFitted pattern center (average of %g patterns):\n',N);
    fprintf(fid,'  xpc = %.4f (X* = %.6f)\n',xpc_ave,xpc_ave_frac);
    fprintf(fid,'  ypc = %.4f (Y* = %.6f)\n',ypc_ave,ypc_ave_frac);
    fprintf(fid,'  L = %.2f (Z* = %.6f)',L_ave,L_ave_frac);
    fprintf(fid,'\n');

    fprintf('Summarized data stored in txt file: %s\n',summarypath);
    fclose(fid);
end

