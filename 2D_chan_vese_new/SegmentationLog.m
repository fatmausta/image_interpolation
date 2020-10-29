function [uu, uu_sq, mv,stime, erriter, num] = SegmentationLog(img, hardIn, hardOut, pdf_in, pdf_out, options)

img_orig = 255 * img;

mask = zeros(size(img));
mask(hardIn) = 1;
%mask(35:end-35,35:end-35) = 1;

ur = mask;
[rows,cols] = size(ur);

step = 0.16;

wi  = 1e2;
wdi = 3;
beta =  15;

num_sq = 1000;
uu_sq = zeros(rows, cols, num_sq+1);

cla;imagesc(img_orig), axis image, axis off; colormap gray; hold on;
contour(ur,[0.5 0.5],'g','LineWidth', 2.5); drawnow;
mv(1) = getframe;
mv(2) = getframe;
uu_sq(:,:,1) = double(ur);
saveas(gcf,'multi_object_rescaled--init.png');

pause(2);

h = fspecial('gaussian',5, 2.0);

img_smooth = imfilter(255.0*img,h);
[dx, dy] =gradient(img_smooth);
grad_mag =  (dx.^2 + dy.^2);

max_grad = max(grad_mag(:));
norm_grad = grad_mag/max_grad;

alpha1 = 0.4; alpha2 = 1.7;  alpha3 = 30;
penalty = alpha1 + alpha2*exp(-alpha3*norm_grad);
penalty_fix = penalty;

bw1  = bwdist(ur,'euclidean');
bw2 = bwdist(1-ur,'euclidean');
bw3 =  bw1 - bw2;
%L = del2(bw3);

[dx, dy] = gradient(bw3);
L2 = (divergence(dx, dy)).^2;

Lmax = max(L2(:));
Lmin = min(L2(:));
L2 = (L2 - Lmin)/(Lmax - Lmin);

penalty = penalty - beta * L2;
penalty(penalty<0) = 1e-3;

IND = floor((255 * img)/options.nbins) + 1;

rf1 =  -log(pdf_in(IND)) ;
rf2 =  -log(pdf_out(IND));

rf1 =   rf1/(max(rf1(:))+eps);
rf2 =   rf2/(max(rf2(:))+eps);

fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .* (rf2 - rf1) * wdi;
fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur).* (rf1 - rf2)* wdi;

fCs(hardIn) = 1e7;
fCt(hardOut) = 1e7;

inpx = double(zeros(rows,cols+1));
inpy = double(zeros(rows+1,cols));
inps = max(fCs, fCt);
inpt = inps;
inpu = double(ur);

varParas = [rows; cols; 300; 1e-5; 0.58; step];

stime = 0;

for i = 1:num_sq
    
    % CPU version
    %  [uu, erriter,num,tt,ps,pt,px,py] = DCMF_mex02(single(penalty), single(fCs), ...
    %      single(fCt), single(varParas), single(inps), single(inpt), single(inpx), ...
    %      single(inpy), single(inpu));
    
    % GPU version
    [uu, erriter,num,tt,ps,pt,px,py] = DCMF_GPU02(single(penalty), single(fCs), ...
        single(fCt), single(varParas), single(inps), single(inpt), single(inpx), ...
        single(inpy), single(inpu));
    stime = tt + stime;
    inpx = px;
    inpy = py;
    inps = ps;
    inpt = pt;
    
    ur = double(uu > 0.5);
    inpu = uu;
    
    cla;colormap(gray);imagesc(img_orig); axis off; axis image; hold on;
    contour(uu,[0.5 0.5],'g','LineWidth', 3); drawnow;
    mv(i+2) = getframe;
    
    bw1  = bwdist(ur,'euclidean');
    bw2 = bwdist(1-ur,'euclidean');
    bw3 =  bw1 - bw2;
    %   L2 = del2(bw3);
    
    [dx, dy] = gradient(bw3);
    L2 = (divergence(dx, dy)).^2;
    Lmax = max(L2(:));
    Lmin = min(L2(:));
    L2 = (L2 - Lmin)/(Lmax - Lmin);
    
    penalty = penalty_fix -  beta * L2;
    penalty(penalty<0) = 1e-2;
    
    
    rf1 =  -log(pdf_in(IND)) ;
    rf2 =  -log(pdf_out(IND));
    
    rf1 =   rf1/(max(rf1(:))+eps);
    rf2 =   rf2/(max(rf2(:))+eps);
    
    fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .* (rf2-rf1) * wdi;
    fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur).* (rf1-rf2)* wdi;
    
    fCs(hardIn) = 1e7;
    fCt(hardOut) = 1e7;
end