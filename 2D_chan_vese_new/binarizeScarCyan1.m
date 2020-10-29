
imOrig = imread('image4--orig.tif');
imcrop = imOrig(175:(175+152),140:(140+164));

figure;imagesc(imcrop); axis image; axis off; colormap gray;

bw = imread('0001.tif');

IMG = size
imgR =imOrig;
imgG = imOrig;
imgB = imOrig;

imgR(bw==128)=250;
imgG(bw==128)=0;
imgB(bw==128) = 0;

imgR(bw==64) = 0;
imgG(bw==64)= 0;
imgB(bw==64) = 250;

IMG(:,:,1) = imgR;
IMG(:,:,2) = imgG;
IMG(:,:,3) = imgB;

figure;imagesc(IMG); axis image; axis off; colormap gray;

