function makebinary_multiimg(patterns, outpath)
%MAKEBINARY_MULTIIMG Take input array of patterns and save as binary .data file
%%% Inputs:
% -patterns: array containing all patterns (height x width x #patterns)
% -outpath: full path to .data file to save as
% 2/21/20


% figure out some things
imagewidth = size(patterns,2);
imageheight = size(patterns,1);
N = size(patterns,3);   % number of patterns


% loop through patterns and write to file
fid = fopen(outpath,'w');       % Open .data file to write   
im = single(patterns);    % Convert data to single precision
im = flipud(im);    % Flip image upside down (TSL has patterns upside down)
       %%% Need above line for EMsoft 4.0, not for EMsoft 4.3
	   
for ii=1:N
    % Make a 1D vector of image
    p = zeros(imageheight*imagewidth,1);     % Initialize vector
    for jj = 1:imageheight
        for kk = 1:imagewidth
            p((jj-1)*imagewidth+kk) = im(jj,kk,ii);
        end
    end
    fwrite(fid,p,'float32');
end


fclose('all');


end

