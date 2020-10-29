function [uu, uu_sq, stime, erriter, num] = SegmentationBh(img, bw_myo, mask,pdf_in, pdf_out, options)
%
%   Function Segmentation_new
%
%   Decription: test the DMCF_mex function for surface evolution

img_orig = img;
%img =  img/max(max(img));

% Creating a new mask
mask = mask;
bw_myo = double(bw_myo)/255;

ur = double(mask)/255;
[rows,cols] = size(ur);

% Parameters for segmentation
wi = 1e3;
wdi = 700;
step = 0.10;

num_sq = 1;
uu_sq = zeros(rows, cols, num_sq+1);

cla;imagesc(img_orig), axis image, axis off; colormap gray; hold on;
contour(ur,[0.5 0.5],'g','LineWidth', 2.5); drawnow;
% mv(1) = getframe;
% mv(2) = getframe;
uu_sq(:,:,1) = double(ur);
%saveas(gcf,'multi_object_rescaled--init.png');

tic;
% Gaussian Filtering
h = fspecial('gaussian',5, 1.5);
img_smooth = imfilter(img,h);
[dx, dy] =gradient(img_smooth);
grad_mag =  (dx.^2 + dy.^2);

max_grad = max(grad_mag(:));
norm_grad = grad_mag/max_grad;

%alpha(1) = 0.3; alpha(2) = 0.1;  alpha(3) = 4;
alpha(1) = 0.27; alpha(2) = 0.1;  alpha(3) = 4;
penalty = alpha(1) + alpha(2)*exp(-alpha(3)*norm_grad);

varParas = [rows; cols; 300; 1e-5; 0.28; step];

% Data term %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samples = 0:1:255; Width =7;

indx = (ur == 1);
Ivals = img (indx);
Ain  =   nnz(indx);
% Intensities scaled to 255
Ivals = Ivals ;
pdf_current_in = ksdensity(Ivals,samples,'Width',Width);

bw = double(bw_myo) - ur;
indx = bw==1;
Aout  =   nnz(indx);
Ivals = img(indx);
Ivals = Ivals ;
pdf_current_out = ksdensity(Ivals,samples,'Width',Width);

IND =  floor(img) + 1;

pdf_in = pdf_in/Ain;
pdf_out = pdf_out/Aout;

p_hat = pdf_current_in/Ain  + eps;
q_hat = pdf_current_out/Aout + eps;

B_in  = (sqrt(pdf_current_in .* pdf_in));
B_out = (sqrt(pdf_current_out .* pdf_out));

B_in  = sum(B_in(:));
B_out = sum(B_out(:));

lookup_in =  (sqrt(pdf_in./p_hat) - B_in) /(2* Ain) ;
lookup_out = (sqrt(pdf_out./q_hat)- B_out)/(2* Aout) ;

rf1 = lookup_in(IND) + B_in/ (2* Ain) ;
rf2 = lookup_out(IND)+ B_out/(2* Aout);

fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .*   (rf1) * wdi;
fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur) .*(rf2 )* wdi;

%fCs =   ur .*   (rf1) * wdi;
%fCt =  (1-ur) .*(rf2 )* wdi;

bw_bg = 1- bw_myo;
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
    
    cla;colormap(gray);imagesc(img_orig); axis off; axis image; hold on;
    contour(uu,[0.5 0.5],'g','LineWidth', 2); drawnow;
    %   mv(i+2) = getframe;
    
    indx = (ur == 1);
    Ivals = img (indx);
    Ain  =   nnz(indx);
    % Intensities scaled to 255
    Ivals = Ivals;
    pdf_current_in = ksdensity(Ivals,samples,'Width',Width);
    
    bw = bw_myo - ur;
    indx = bw==1;
    Aout  =   nnz(indx);
    Ivals = img(indx);
    Ivals = Ivals ;
    pdf_current_out = ksdensity(Ivals,samples,'Width',Width);
    %
    IND =  floor( img) + 1;
    %
    pdf_in = pdf_in/Ain;
    pdf_out = pdf_out/Aout;
    %
    p_hat = pdf_current_in/Ain  + eps;
    q_hat = pdf_current_out/Aout + eps;
    %
    B_in  = (sqrt(pdf_current_in .* pdf_in));
    B_out = (sqrt(pdf_current_out .* pdf_out));
    %
    B_in  = sum(B_in(:));
    B_out = sum(B_out(:));
    %
    lookup_in =  (sqrt(pdf_in./p_hat) - B_in) /(2* Ain) ;
    lookup_out = (sqrt(pdf_out./q_hat)- B_out)/(2* Aout) ;
    %
    rf1 = lookup_in(IND) + B_in/ (2* Ain) ;
    rf2 = lookup_out(IND)+ B_out/(2* Aout);
    %
    fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .*   (rf1) * wdi;
    fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur) .*(rf2 )* wdi;
    %
    fCt(logical(bw_bg)) = 1e7;
end