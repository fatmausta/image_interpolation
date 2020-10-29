function [uu, uu_sq, mv,stime, erriter, num] = SegmentationCV(img, bw_myo, mask)
%
%   Function Segmentation_new
%
%   Decription: test the DMCF_mex function for surface evolution

options.GPU = 0;

img_orig = img;
img = img/max(max(img));

% Creating a new mask
mask = mask;

ur = double(mask)/255;
[rows,cols] = size(ur);

% Parameters for segmentation
wi = 1e6;
wdi = 300;
step = 0.10;

num_sq = 1;
uu_sq = zeros(rows, cols, num_sq+1);

cla;imagesc(img_orig), axis image, axis off; colormap gray; hold on;
contour(ur,[0.5 0.5],'g','LineWidth', 2.5); drawnow;
mv(1) = getframe;
mv(2) = getframe;
uu_sq(:,:,1) = double(ur);
%saveas(gcf,'multi_object_rescaled--init.png');

tic;
% Gaussian Filtering
%h = fspecial('gaussian',7, 1.0);
%img_smooth = imfilter(255.0*img,h);
[dx, dy] =gradient(img);
grad_mag =  (dx.^2 + dy.^2);

max_grad = max(grad_mag(:));
norm_grad = grad_mag/max_grad;

alpha(1) = 0.1; alpha(2) = 0.2;  alpha(3) = 2;
penalty = alpha(1) + alpha(2)*exp(-alpha(3)*norm_grad);

varParas = [rows; cols; 300; 1e-5; 0.38; step];

% Inside
indx = find(ur == 1);
Ain  =   size(indx,1);
Ivals = img (indx);
% Intensities scaled to 255;
Iin = mean(Ivals);

% Outside
indx = find(ur == 0);
% Volume computation
Aout  =   size(indx,1);
Ivals = img(indx);
% Intensities scaled to 255
Iout = mean(Ivals);

rf1 = (img - Iin).^2 ;
rf2 = (img - Iout).^2;

rf1 =   rf1/(max(rf1(:))+1e-7);
rf2 =   rf2/(max(rf2(:))+1e-7);

fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .* (rf2-rf1) * wdi;
fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur).* (rf1-rf2)* wdi;

% fCs = ur .* rf1 * wdi;
% fCt = (1-ur) .*rf2* wdi;

bw_bg = 1 - (double(bw_myo)/255);
fCt(logical(bw_bg)) = 1e7;

inpx = double(zeros(rows,cols+1));
inpy = double(zeros(rows+1,cols));
inps = max(fCs, fCt);
inpt = inps;
inpu = double(ur);

stime = 0;

for i = 1:num_sq
    
    if(options.GPU == 0)
        % CPU version
        [uu, erriter,num,tt,ps,pt,px,py] = DCMF_mex02(single(penalty), single(fCs), ...
            single(fCt), single(varParas), single(inps), single(inpt), single(inpx), ...
            single(inpy), single(inpu));
    elseif(options.GPU == 1)
        % GPU version
        [uu, erriter,num,tt,ps,pt,px,py] = DCMF_GPU02(single(penalty), single(fCs), ...
            single(fCt), single(varParas), single(inps), single(inpt), single(inpx), ...
            single(inpy), single(inpu));
    end
    
    stime = tt + stime;
    inpx = px;
    inpy = py;
    inps = ps;
    inpt = pt;
    
    ur = double(uu > 0.5);
    inpu = uu;
    
    organ = 'multi';
    cla;colormap(gray);imagesc(img_orig); axis off; axis image; hold on;
    contour(uu,[0.5 0.5],'g','LineWidth', 2); drawnow;
    mv(i+2) = getframe;
%     imwrite(mv(i+1).cdata,[organ '\' num2str(i+1) '.png' ]);
    
    %  elapsedTime = toc;
    %  totalTime = totalTime+elapsedTime;
    %  disp(['Total time: ' num2str(totalTime)]);
    %  key = input('Continue or not: ', 's');
    %   if isempty(key)
    %  else
    %   break;
    %   end
     
    % Inside
    indx = find(ur == 1);
    Ain  =   size(indx,1);
    Ivals = img (indx);
    % Intensities scaled to 255
    Iin = mean(Ivals);
        
    % Outside
    indx = find(ur == 0);
    % Volume computation
    Aout  =   size(indx,1);
    Ivals = img(indx);
    % Intensities scaled to 255
    Iout = mean(Ivals);
    
    rf1 = (img - Iin).^2 ;
    rf2 = (img - Iout).^2;
    
    rf1 =   rf1/(max(rf1(:))+1e-7);
    rf2 =   rf2/(max(rf2(:))+1e-7);
    
    fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .* (rf2-rf1) * wdi;
    fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur).* (rf1-rf2)* wdi;
    
    fCt(logical(bw_bg)) = 1e7;
end