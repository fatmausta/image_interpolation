% Main function for Interpolating scar tissue
% Author: Eranga Ukwatta
% Date: feb 1, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting the figure location
set(0,'defaultfigureposition',[400 400 600 600]);

% Choices
isd = 8; isd_mm = 1/isd;

interpMethod = 'spline';
LogOddsSpace = 'gaussian';

load('DE');

% 3 for CMF with Bh  
% 5 FWHM 
% 7 Region Growing  
% 8 for STRM+1SD
% 9 for STRM+2SD
% 10 for STRM+3SD
% 11 for STRM+4SD

method = 7;

INo = [26 29 54];

for i = 2:12%1:3
    
    %i = INo(i);
    slices = DE{i};
    imgSparse = [];
    jj = 0;
    
    for j = 1:size(slices,2)
        bw = imread(['results/' num2str(i) '/m'  num2str(method) '-'  num2str(j)  '-bw.png']);
        imgSparse(:,:,j) = bw;
    end
    
    % Apply Log Odds based interpolation technique
    imageOut = interpBinImageLogOdds(imgSparse,1:size(imgSparse,2),1:size(imgSparse,1), 1:isd_mm:size(imgSparse,3),0,interpMethod,LogOddsSpace);
    imageOut = logistic_fn(imageOut)>0.5;
    
    figure; clf; isosurface(imageOut,0.5); axis image; axis off; drawnow;
    
    %metaImageWrite(uint8(255*imageOut),['results/'  num2str(i) '/m'  num2str(method)  '-bw_a.mhd']);
    
end
