function makebinary_singleimg(img_exp_path, outpath, data)
%MAKEBINARY_SINGLEIMG Load image and save as binary .data file
%%% Inputs:
% -img_exp_path: full path to image file to read in
% -outpath: full path to .data file to save as
% -data: struct containing binning, numsx, numsy
%
% Original: 12/16/19 (Edward Pang, MIT)


% extract data
binning = data.binning;
numsx = data.numsx;
numsy = data.numsy;

% load experimental pattern
img_exp = imread(img_exp_path);
if size(img_exp,3)>1
    img_exp = img_exp(:,:,1);   % remove extra channels if necessary
end

% figure out image sizes
if size(img_exp,2) ~= numsx/binning || size(img_exp,1) ~= numsy/binning
	error('Image size does not match numsy/binning x numsx/binning.');
end
imagewidth = numsx/binning;        % Width (px) of patterns
imageheight = numsy/binning;        % Height (px) of patterns


% Read in image
fid = fopen(outpath,'w');       % Open .data file to write   
im = single(img_exp);    % Convert data to single precision
im = flipud(im);    % Flip image upside down (TSL has patterns upside down)
	%%% Need above line for EMsoft 4.0, not for EMsoft 4.3
% Make a 1D vector of image
p = zeros(imageheight*imagewidth,1);     % Initialize vector
for jj = 1:imageheight
    for kk = 1:imagewidth
        p((jj-1)*imagewidth+kk) = im(jj,kk);
    end
end
fwrite(fid,p,'float32');
fclose('all');


end

