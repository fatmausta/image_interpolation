function [uu, uu_sq, stime, erriter, num] = SegmentationLog2(img, bw_myo, mask, pdf_in, pdf_out)
%
%   Function Segmentation_new
%
%   Decription: test the DMCF_mex function for surface evolution

img_orig = img;
img =  img/max(max(img));

% Creating a new mask
mask = mask;

ur = double(mask)/255;
[rows,cols] = size(ur);

% Parameters for segmentation
%wi = 1e6;
wi = 1e4;
wdi = 3;
step = 0.1;

num_sq = 1;
uu_sq = zeros(rows, cols, num_sq+1);

%cla;imagesc(img_orig), axis image, axis off; colormap gray; hold on;
%contour(ur,[0.5 0.5],'g','LineWidth', 2.5); drawnow;
%mv(1) = getframe;
%mv(2) = getframe;
uu_sq(:,:,1) = double(ur);
%saveas(gcf,'multi_object_rescaled--init.png');

tic;
% Gaussian Filtering
h = fspecial('gaussian',5, 1.0);
img_smooth = imfilter(255.0*img,h);
[dx, dy] =gradient(img_smooth);
grad_mag =  (dx.^2 + dy.^2);

max_grad = max(grad_mag(:));
norm_grad = grad_mag/max_grad;

alpha(1) = 0.05; alpha(2) = 0.1;  alpha(3) = 3;
penalty = alpha(1) + alpha(2)*exp(-alpha(3)*norm_grad);

varParas = [rows; cols; 300; 1e-5; 0.38; step];

IND = floor((255 * img)) + 1;

rf1 =  abs(-log(pdf_in(IND))) ;
rf2 =  abs(-log(pdf_out(IND)));

%rf1 =   rf1/(max(rf1(:))+eps);
%rf2 =   rf2/(max(rf2(:))+eps);

%fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .* (rf2-rf1) * wdi;
%fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur).* (rf1-rf2)* wdi;

 fCs = (rf2) * wdi;
 fCt = (rf1)* wdi;

bw_bg = 1- (double(bw_myo)/255);
fCt(logical(bw_bg)) = 1e7;

inpx = double(zeros(rows,cols+1));
inpy = double(zeros(rows+1,cols));
inps = max(fCs, fCt);
inpt = inps;
inpu = double(ur);

stime = 0;

for i = 1:num_sq
    
    % CPU version
      [uu, erriter,num,tt,ps,pt,px,py] = DCMF_mex02(single(penalty), single(fCs), ...
          single(fCt), single(varParas), single(inps), single(inpt), single(inpx), ...
          single(inpy), single(inpu));
    
    % GPU version
%      [uu, erriter,num,tt,ps,pt,px,py] = DCMF_GPU02(single(penalty), single(fCs), ...
%        single(fCt), single(varParas), single(inps), single(inpt), single(inpx), ...
%        single(inpy), single(inpu));

    stime = tt + stime;
    inpx = px;
    inpy = py;
    inps = ps;
    inpt = pt;
    
    ur = double(uu >= 0.5);
    inpu = uu;
    
    cla;colormap(gray);imagesc(img_orig); axis off; axis image; hold on;
    contour(uu,[0.5 0.5],'g','LineWidth', 2); drawnow;
 %   mv(i+2) = getframe;
    
    
    rf1 =  -log(pdf_in(IND)) ;
    rf2 =  -log(pdf_out(IND));
    
    rf1 =   rf1/(max(rf1(:))+eps);
    rf2 =   rf2/(max(rf2(:))+eps);
    
    fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .* (rf2-rf1) * wdi;
    fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur).* (rf1-rf2)* wdi;
    
    fCt(logical(bw_bg)) = 1e7;
end