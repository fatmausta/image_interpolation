function distOut = interpBinImageLogOdds(imageIn, xVals, yVals, zVals, smoothingIterations, interpMethod, mapType)
%
% This function interpolates a given 3D binary image using distance
% transform to preserve shape
%
% Input:
% imageIn       -   The input binary image
% inRes         -   Resolution of the input image: [resX resY resZ]
% xVals         -   The x values for the voxels of the new image. The
%                   original image is assumed to have x values 1 to
%                   size(image, 2)
%

% The interpolation parameters
if(~exist(interpMethod))
    interpMethod = 'linear';
end
if (~exist('smoothingIterations'))
    smoothingIterations = 0;
end
if (~exist('Thres'))
    Thres = 0;
end
if(~exist('mapType'))
    mapType = 'gaussian';
end

% Ensure binary
imageIn(imageIn>0) = 1;
imageIn(imageIn<0) = 0;

% Image dimensions
dimZin = size(imageIn, 3);
distIn = zeros(size(imageIn));

if(strcmp(mapType,'gaussian'))
    % Smoothing using a Gaussian function
   % h = fspecial('gaussian',9,3);
    h = fspecial('gaussian',9,4);
     
  %  distIn = imgaussian(double(imageIn),1.5,7);
    
    for i = 1:dimZin
        distIn(:,:,i) = imfilter(imageIn(:,:,i),h);
        temp = distIn(:,:,i);
        temp = temp/(max(temp(:)+1e-9));  % 1e-13
        %distIn(:,:,i) = temp;
        temp = log((temp)./(1-temp));
        %    temp = 1./(1+ exp(-temp));
        temp(temp==-Inf) =   -10;
        distIn(:,:,i) = temp;
    end
elseif(strcmp(mapType,'signedDist'))
    % Calculate the signed distance transform of each of the slices
    disp('Extracting perimeter and computing distance transform ...');
    
    for i = 1:dimZin
        if (nnz(imageIn(:, :, i)) == 0)
            distIn(:, :, i) = numel(zVals)/abs(zVals(end)-zVals(1));
            continue;
        end
        perimSlice = bwperim(imageIn(:, :, i));
        distSlice = bwdist(perimSlice);
        distSlice(imageIn(:, :, i)>0) = -distSlice(imageIn(:, :, i)>0);
        distIn(:, :, i) = distSlice;
        
        %distSlice = bwdist(imageIn(:,:,i));
        %distSlice1= bwdist(1-imageIn(:,:,i));
        %distSlice = -distSlice + distSlice1;
        %distIn(:, :, i) = distSlice;
    end
end

% Interpolate the distance transform image
disp('Interpolating ...');
[Xi, Yi, Zi] = meshgrid(xVals, yVals, zVals);
distOut= interp3(distIn, Xi, Yi, Zi, interpMethod);

% Smooth the distance image
disp('Smoothing the distance image ...');
for i = 1:smoothingIterations
    distOut = smooth3(distOut);
end

% Threshold the dist image
%imageOut = (distOut <= Thres);
%imageOut = distIn;