function loadpattern(data)
%LOADPATTERN Load pattern image file, bin, and save as binary .data file
%%% Inputs:
% -data: struct containing location of input/output files and other relevant info
% Assumes image is as captured from a TSL/EDAX system. Pattern recorded on a Bruker 
% or Oxford system may require a flip/rotation (not currently implemented, but you 
% can manually do this to agree with TSL coordinate system)
%
% Original: 7/24/19 (Edward Pang, MIT)
% Change log:
% -10/10/19 ELP: Comment out flipup to make coordinate system agree with EMsoft 4.3
% -2/20/20 ELP: throw error if image size not equal to specified nums/binning


    % Extract parameters from struct 'data'
    homepath = data.homepath;
    path = data.path;
    tmppath = data.tmppath;
    imgprebinned = data.imgprebinned;
    numsx = data.numsx;
    numsy = data.numsy;
    binning = data.binning;
    img_exp_name = data.img_exp_name;
    outname = data.img_exp;
    

    % Figure out paths
    img_exp_path = fullfile(homepath,path,img_exp_name);
    outpath = fullfile(homepath,tmppath,outname);

    
    % load experimental pattern
    img_exp = imread(img_exp_path);
    if size(img_exp,3)>1
        img_exp = img_exp(:,:,1);   % remove extra channels if necessary
    end

    % bin it down
    img_exp_bin = zeros(numsy/binning,numsx/binning,'uint8');   % initialize
    if imgprebinned==0
        for kk=0:binning:numsy-binning
            for ll=0:binning:numsx-binning
                img_exp_bin(kk/binning+1,ll/binning+1) = sum(sum(img_exp(kk+1:kk+binning,ll+1:ll+binning)))/binning^2;
            end
        end
    else
        img_exp_bin = img_exp;
    end


    % figure out some parameters for loading in images
    imagewidth = numsx/binning;        % Width (px) of patterns
    imageheight = numsy/binning;        % Height (px) of patterns
    
    if imagewidth~=size(img_exp_bin,2) || imageheight~=size(img_exp_bin,1)
        error('Image size does not match nums/binning.');
    end

    % Read in image
    fid = fopen(outpath,'w');       % Open .data file to write   
    im = single(img_exp_bin);    % Convert data to single precision
    %im = flipud(im);    % Flip image upside down (TSL has patterns upside down)
        %%% Need above line for EMsoft 4.0, not for EMsoft 4.3
    % Make a 1D vector of image
    P = zeros(imageheight*imagewidth,1);     % Initialize vector
    for jj = 1:imageheight
        for kk = 1:imagewidth
            P((jj-1)*imagewidth+kk) = im(jj,kk);
        end
    end
    fwrite(fid,P,'float32');
    fclose('all');

end

