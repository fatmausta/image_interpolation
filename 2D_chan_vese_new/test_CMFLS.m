function [uu, uu_sq, mv, erriter, num] = test_CMFLS
%
%   Function test_CMFLS
%
%   Decription: test the DMCF_mex function for surface evolution
%
%   Author: Jing Yuan, UWO, London
%
ur = zeros(300,290);
%ur=double(ur)/255;
ur(100:200,100:150) = 1;
ur(150:200,100:200) = 1;
[rows,cols] = size(ur);

wi = 60;
num_sq = 30;
uu_sq = zeros(rows, cols, num_sq+1);

imagesc(ur), axis image, axis off;
mv(1) = getframe;
uu_sq(:,:,1) = double(ur);

penalty = 0.5*ones(rows,cols);
varParas = [rows; cols; 300; 1e-5; 0.43; 0.16];

fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi;
fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi;
inpx = double(zeros(rows,cols+1));
inpy = double(zeros(rows+1,cols));
inps = max(fCs, fCt);
inpt = inps;
inpu = double(ur);

% CPU version
%[uu, erriter,num,tt, ps, pt, px, py] = DCMF_mex02(single(penalty), single(fCs), single(fCt), ... 
%         single(varParas), single(inps), single(inpt), single(inpx), single(inpy), single(inpu));
     
% GPU version
 [uu, erriter,num,tt, ps, pt, px, py] = DCMF_GPU02(single(penalty), single(fCs), single(fCt), ... 
          single(varParas), single(inps), single(inpt), single(inpx), single(inpy), single(inpu));
     

uu_sq(:,:,2) = double(uu>0.5);
imagesc(uu_sq(:,:,2)), axis image, axis off;
mv(2) = getframe;

for i = 1:num_sq-1

     inpx = px;
     inpy = py;
     inps = ps;
     inpt = pt;

     ur = double(uu > 0.5);
     inpu = uu;
    
     fCs = (ur.*(bwdist(1-ur, 'euclidean')))/wi;
     fCt = ((1-ur).*(bwdist(ur,'euclidean')))/wi;
    
     % CPU version   
%     [uu, erriter,num,tt,ps,pt,px,py] = DCMF_mex02(single(penalty), single(fCs), ...
%        single(fCt), single(varParas), single(inps), single(inpt), single(inpx), ...
%        single(inpy), single(inpu));
    
    % GPU version
     [uu, erriter,num,tt,ps,pt,px,py] = DCMF_GPU02(single(penalty), single(fCs), ...
         single(fCt), single(varParas), single(inps), single(inpt), single(inpx), ...
         single(inpy), single(inpu));

    uu_sq(:,:,i+2) = double(uu>0.5);
    imagesc(uu_sq(:,:,i+2)), axis image, axis off;
    mv(i+2) = getframe;
    
end
