addpath('imMinkowski');

% 3 for CMF with Bh
% 5 FWHM
% 7 Region Growing %
% 8 for STRM+1SD
% 9 for STRM+2SD
% 10 for STRM+3SD
% 11 for STRM+4SD

method = 5;

load('Res');

for i = 1:60
    
    bw_a =     metaImageRead(['results/'  num2str(i) '/m'  num2str(method)  '-bw_a.mhd']);
    bw_m =     metaImageRead(['results/'  num2str(i) '/'  '-bw_m.mhd']);
    
    bw_a = logical(bw_a);
    bw_m = logical(bw_m);
    
    Dice(i) = Dice_Index(logical(bw_a), (bw_m));
    
    % Root mean square error
    rmse = RMSE3D(bw_a,bw_m, Res(i));
    Dist(i) = mean(rmse);
    
    [f1, autoSet] = isosurface(bw_a,0);
    [f2, refSet]   = isosurface(bw_m,0);
    
    indx = sub2ind(size(bw_a),autoSet(:,2),autoSet(:,1),autoSet(:,3));
    imgSa = zeros(size(bw_a)); imgSa(indx) = 1;
    
    indx      = sub2ind(size(bw_m),refSet(:,2),refSet(:,1),refSet(:,3));
    imgSm = zeros(size(bw_m)); imgSm(indx) = 1;
    
    area_m(i) =  imSurface(logical(imgSm),[1, 1, 1]);
    area_a(i)  =  imSurface(logical(imgSa), [1, 1, 1]);
    
    [chi_m(i), labels_m] = imEuler3d(imgSm,26);
    [chi_a(i), labels_a] = imEuler3d(imgSa,26);
    
    
    options.curvature_smoothing = 0;
    options.verb = 0;
 %   [Umin,Umax,Cmin,Cmax_A,Cmean_A,Cgauss_A,Normal]           = compute_curvature(autoSet,f1,options);
    
    %[Umin_M,Umax_M,Cmin_M,Cmax_M,Cmean_M,Cgauss_M,Normal_M] = compute_curvature(refSet,f2,options);
    
    %BA = hist(Cmax_A,-0.6:0.02:0.6);
    %BB = hist(Cmax_M,-0.6:0.02:0.6);
    
    %   options.curvature_smoothing =5;
   % BA = BA/sum(BA);
  %  BB = BB/sum(BB);
    
 %   Bh(i) = sum(sqrt(BA.*BB))
end


save(['Dice-m' num2str(method)],'Dice');
save(['Dist-m' num2str(method)],'Dist');

save(['area-m' ],'area_m');
save(['area-m' num2str(method)],'area_a');

save(['chi-m' ],'chi_m');
save(['chi-m' num2str(method)],'chi_a');


%save(['Bh-m'  num2str(method)],'Bh');
