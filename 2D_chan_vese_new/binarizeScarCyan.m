
%
img = imread('Scar.png');
figure;imagesc(img); axis image; axis off; colormap gray;



indCyan  = find(img(:,:,1) ==0  & img(:,:,2) ==255 & img(:,:,3) ==255);

[m, n, p] = size(img);
Contour = zeros(m, n);
Contour(indCyan) = 1;

[B,L,N,A] = bwboundaries(Contour, 8, 'holes');
[r,~] = find(A(:,N+1:end));
[rr,~] = find(A(:,r));
idx = setdiff(1:numel(B), [r(:);rr(:)]);
bw = ismember(L,idx);

imwrite(bw,'bwCyan.tif');



I = img(:,:,2);


I = I .*uint8(bw);

cp = find(I >60);

imgR = img(:,:,1);
imgG = img(:,:,2);
imgB = img(:,:,3);

imgR(cp)=250;
imgG(cp)=0;
imgB(cp) = 0;

img = zeros(size(img));
IMG(:,:,1) = imgR;
IMG(:,:,2) = imgG;
IMG(:,:,3) = imgB;

figure;imagesc(IMG); axis image; axis off; colormap gray;

% Border zone
bw(cp) = 0;
bp = find(bw==1);

imgR(bp) = 0;
imgG(bp)= 0;
imgB(bp) = 250;

IMG(:,:,1) = imgR;
IMG(:,:,2) = imgG;
IMG(:,:,3) = imgB;

figure;imagesc(IMG); axis image; axis off; colormap gray;

