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
INo = [];

for   i =  1:1
    
    slices = DE{i};
    imgSparse = [];
    jj = 0;
    
    for j = 1:size(slices,2)
        bw = imread(['results/' num2str(i) '/man-'  num2str(j)  '-bw.png']);
        imgSparse(:,:,j) = bw;
    end
    
    % Apply Log Odds based interpolation technique
    imageOut = interpBinImageLogOdds(imgSparse,1:size(imgSparse,2),1:size(imgSparse,1), 1:isd_mm:size(imgSparse,3),0,interpMethod,LogOddsSpace);
    imageOut = logistic_fn(imageOut)>0.5;
    
    clf;isosurface(imageOut,0.5);axis image;axis off;drawnow;
    
    metaImageWrite(uint8(255*imageOut),['results/'  num2str(i) '/'    '-bw_m'   '.mhd']);
end
