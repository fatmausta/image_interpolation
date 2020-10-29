function [uu, uu_sq, mv, erriter, num] = Segmentation(img,mask,pdf_prior_in,pdf_prior_out,forePixels,backPixels, IM)
%
%   Function Segmentation_new
%
%   Decription: test the DMCF_mex function for surface evolution
%
%   Author:   Jing Yuan, UWO, London
%   Modified: Eranga Ukwatta

img_orig = img;
img =  img/max(max(img));

% Creating a new mask
mask = zeros(size(img));
mask(25:end-25,25:end-25) = 1;

ur = mask;
[rows,cols] = size(ur);

% Parameters for segmentation
wi = 1e6;
wdi = 300;

% Baby brains
%wdi = 2000;

samples = 0:1:255;
Width = 10;
step = 0.14;

%Prostate
%organ = 'prostate';
%wi = 1e5;
%wdi = 20;

% Multi-object
wi = 1e4;
wdi = 5200;

% new Liver 
%wi = 1e6;
%wdi = 3500;

% Zebra
wi = 5e3;
wdi = 11500;

%Liver
%organ = 'liver';
%wdi = 500;

% bacteria
%wi = 10000;
%wdi = 12000;
%step = 0.08;

totalTime = 0;

num_sq = 38;
uu_sq = zeros(rows, cols, num_sq+1);

cla;imagesc(img_orig), axis image, axis off; colormap gray;
mv(1) = getframe;
mv(2) = getframe;
uu_sq(:,:,1) = double(ur);

tic;
% Gaussian Filtering

% Zebra
h = fspecial('gaussian',5, 0.1);
% multi-object
%h = fspecial('gaussian',5, 0.4);
img_smooth = imfilter(255.0*img,h);
[dx, dy] =gradient(img_smooth);
grad_mag =  (dx.^2 + dy.^2);

max_grad = max(grad_mag(:));
norm_grad = grad_mag/max_grad;
%img = img_smooth./(max(max(img_smooth)));

%penalty = 0.3*ones(rows,cols); 
% Liver 2
%alpha1 = 0.08; alpha2 = 0.4;  alpha3 = 40;
% Zebra
alpha1 = 0.2; alpha2 = 0.7;  alpha3 = 23;

% multiobject
%alpha1 = 0.08; alpha2 = 0.3;  alpha3 = 15;

penalty = alpha1 + alpha2*exp(-alpha3*norm_grad);

varParas = [rows; cols; 300; 1e-5; 0.58; step];
p     = pdf_prior_in;
q     = pdf_prior_out;

% Inside
indx = find(ur == 1);
Ain  =   size(indx,1);
Ivals = img (indx);
% Intensities scaled to 255
Ivals = Ivals * 255;
pdf_current_in = ksdensity(Ivals,samples,'Width',Width);
%figure,plot(samples,pdf_current_in);axis([0 150 0 0.05]); hold on; grid on;

% Outside
indx = find(ur == 0);
% Volume computation
Aout  =   size(indx,1);
Ivals = img(indx);
% Intensities scaled to 255
Ivals = Ivals * 255;
pdf_current_out = ksdensity(Ivals,samples,'Width',Width);
%plot(samples,pdf_current_out);axis([0 150 0 0.05]);

p_hat = pdf_current_in;
q_hat = pdf_current_out;

IND =  floor(255* img) + 1;

p_hat =  p_hat + eps;
q_hat =  q_hat + eps;
B_in =  sum(sqrt(pdf_current_in.*p));
B_out = sum(sqrt(pdf_current_out.*q));
lookup_in =    (sqrt(p./p_hat) - B_in )/(2*Ain) ;
lookup_out =   (sqrt(q./q_hat) - B_out)/(2*Aout);

rf1 = lookup_in(IND)  ;
rf2 = lookup_out(IND) ;

fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .* (rf1) * wdi;
fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur) .* (rf2)* wdi;

% fCs = ur .* rf1 * wdi;
% fCt = (1-ur) .*rf2* wdi;

fCs(forePixels) = 1e7;
fCt(backPixels) = 1e7;

inpx = double(zeros(rows,cols+1));
inpy = double(zeros(rows+1,cols));
inps = max(fCs, fCt);
inpt = inps;
inpu = double(ur);

for i = 1:num_sq
    
    % CPU version
  %   [uu, erriter,num,tt,ps,pt,px,py] = DCMF_mex02(single(penalty), single(fCs), ...
  %       single(fCt), single(varParas), single(inps), single(inpt), single(inpx), ...
  %       single(inpy), single(inpu));
    
    % GPU version
   [uu, erriter,num,tt,ps,pt,px,py] = DCMF_GPU02(single(penalty), single(fCs), ...
        single(fCt), single(varParas), single(inps), single(inpt), single(inpx), ...
        single(inpy), single(inpu));
    
    inpx = px;
    inpy = py;
    inps = ps;
    inpt = pt;
    
    ur = double(uu > 0.5);
    inpu = uu;
    
    cla;colormap(gray);imagesc(img_orig); axis off; axis image; hold on;
    contour(uu,[0.5 0.5],'g','LineWidth', 2); drawnow;
    mv(i+2) = getframe;
%s    imwrite(mv(i+1).cdata,[organ '\' num2str(i+1) '.png' ]);
    
    elapsedTime = toc;
    totalTime = totalTime+elapsedTime;
    disp(['Total time: ' num2str(totalTime)]);
    %  key = input('Continue or not: ', 's');
    %   if isempty(key)
    %  else
    %   break;
    %   end
    
    tic;
    
    % Inside
    indx = find(ur == 1);
    Ain  =   size(indx,1);
    Ivals = img (indx);
    % Intensities scaled to 255
    Ivals = Ivals * 255.0;
    pdf_current_in = ksdensity(Ivals,samples,'Width',Width);
    %figure,plot(samples,pdf_current_in);axis([0 150 0 0.05]); hold on; grid on;
    
    % Outside
    indx = find(ur == 0);
    % Volume computation
    Aout  =   size(indx,1);
    Ivals = img(indx);
    % Intensities scaled to 255
    Ivals = Ivals * 255.0;
    pdf_current_out = ksdensity(Ivals,samples,'Width',Width);
    %plot(samples,pdf_current_out);axis([0 150 0 0.05]);
    
    p_hat = pdf_current_in;
    q_hat = pdf_current_out;
    
    IND =  floor(255.0* img) + 1;
    
    p_hat =  p_hat + eps;
    q_hat =  q_hat + eps;
    B_in =  sum(sqrt(pdf_current_in.*p))/Ain;
    B_out = sum(sqrt(pdf_current_out.*q))/Aout;
    
    lookup_in =    (sqrt(p./p_hat) - B_in )/(2*Ain) ;
    lookup_out =   (sqrt(q./q_hat)   - B_out)/(2*Aout);
    
    rf1 = lookup_in(IND)  ;
    rf2 = lookup_out(IND) ;
    
    fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .* ( rf1) * wdi;
    fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur) .*( rf2)* wdi;
    
    fCs(forePixels) = 1e7;
    fCt(backPixels) = 1e7;
end