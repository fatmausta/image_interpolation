function [uu, uu_sq, mv, erriter, num] = Segmentation(img,mask,pdf_prior_in,pdf_prior_out,forePixels,backPixels)
%
%   Function Segmentation_new
%
%   Decription: test the DMCF_mex function for surface evolution
%
%   Author:   Jing Yuan, UWO, London
%   Modified: Eranga Ukwatta

img =  img/max(max(img));

% Creating a new mask
%mask = zeros(size(img));
%mask(15:end-15,15:end-15) = 1;

ur = mask;
[rows,cols] = size(ur);

% Parameters for segmentation
bh = 1;
wi = 1e5;
wdi = 200;
alpha = 1;
Beta = 1;
samples = 0:1:255;
Width = 3;

num_sq = 75;
uu_sq = zeros(rows, cols, num_sq+1);

imagesc(img), axis image, axis off;
mv(1) = getframe;
uu_sq(:,:,1) = double(ur);

% Gaussian Filtering
h = fspecial('gaussian',5, 1.0);
img_smooth = imfilter(255.0*img,h);
[dx, dy] =gradient(img_smooth);
grad_mag =  sqrt(dx.^2 + dy.^2);
img = img_smooth./(max(max(img_smooth)));

penalty = 0.3*ones(rows,cols);
%penalty =  alpha ./(1 + Beta * grad_mag);

varParas = [rows; cols; 300; 1e-5; 0.38; 0.16];
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

IND =  round(255* img) + 1;

if(bh == 1)   % Bhattacharyya energy
    p_hat =  p_hat + eps;
    q_hat =  q_hat + eps;
    B_in =  sum(sqrt(pdf_current_in.*p));
    B_out = sum(sqrt(pdf_current_out.*q));
    lookup_in =    (sqrt(p./p_hat) )/(2*Ain) ;
    lookup_out =   (sqrt(q./q_hat) )/(2*Aout);
else         % KL divergence
    p_hat =  p_hat + eps;
    q_hat =  q_hat + eps;
    %lookup_in =  (1 + log(p_hat./p) - p./p_hat)/(Ain)   ;
    %lookup_out = (1 + log(q_hat./q) - q./q_hat)/( Aout) ;
    lookup_in =  (1 + log(p./p_hat)) /(Ain)   ;
    lookup_out = (1 + log(q./q_hat)) /( Aout) ;
end

rf1 = lookup_in(IND) - B_in/(2 * Ain) ;
rf2 = lookup_out(IND)- B_out/(2* Aout) ;

fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .* (rf1) * wdi;
fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur) .* (rf2)* wdi;

% fCs = ur .* rf1 * wdi;
% fCt = (1-ur) .*rf2* wdi;

fCs(forePixels) = 1e8;
fCt(backPixels) = 1e8;

inpx = double(zeros(rows,cols+1));
inpy = double(zeros(rows+1,cols));
inps = max(fCs, fCt);
inpt = inps;
inpu = double(ur);

for i = 1:num_sq
    
    % CPU version
    % [uu, erriter,num,tt,ps,pt,px,py] = DCMF_mex02(single(penalty), single(fCs), ...
    %     single(fCt), single(varParas), single(inps), single(inpt), single(inpx), ...
    %     single(inpy), single(inpu));
    
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
    
    clf;colormap(gray);imagesc(img); axis off; axis image; hold on;
    contour(uu,[0.5 0.5],'g','LineWidth', 2); drawnow;
    mv(i+1) = getframe;
    
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
    
    IND =  round(255* img) + 1;
    
    if(bh == 1)   % Bhattacharyya energy
        p_hat =  p_hat + eps;
        q_hat =  q_hat + eps;
        B_in =  sum(sqrt(pdf_current_in.*p));
        B_out = sum(sqrt(pdf_current_out.*q));
       lookup_in =    (sqrt(p./p_hat)  )/(2*Ain) ;
    lookup_out =   (sqrt(q./q_hat)   )/(2*Aout);
    else         % KL divergence
        p_hat =  p_hat + eps;
        q_hat =  q_hat + eps;
        %  lookup_in =  (1 + log(p_hat./p) - p./p_hat)/Ain   ;
        %  lookup_out = (1 + log(q_hat./q) - q./q_hat)/Aout ;
        lookup_in =  (1 + log(p./p_hat) )  /(Ain)   ;
        lookup_out = (1 + log(q./q_hat) ) /( Aout) ;
    end
    
    
  rf1 = lookup_in(IND) - B_in/(2 * Ain) ;
rf2 = lookup_out(IND)- B_out/(2* Aout) ;
    
    fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi  + ur .* ( rf1) * wdi;
    fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi + (1-ur) .*(rf2)* wdi;
    
    %     fCs = ur .* rf1 * wdi;
    %     fCt = (1-ur) .*rf2* wdi;
    
    fCs(forePixels) =  1e8;
    fCt(backPixels) = 1e8;
    
end
