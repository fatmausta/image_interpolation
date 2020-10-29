addpath('imMinkowski');

% 3 for CMF with Bh
% 5 FWHM
% 7 Region Growing %
% 8 for STRM+1SD
% 9 for STRM+2SD
% 10 for STRM+3SD
% 11 for STRM+4SD

method = 5;

for i = 2:56
    
    bw_a =     metaImageRead(['results/'  num2str(i) '/m'  num2str(method)  '-bw_a.mhd']);
    bw_m =     metaImageRead(['results/'  num2str(i) '/'  '-bw_m.mhd']);
    
  %  bw_a = logical(bw_a);
  %  bw_m = logical(bw_m);
    
    [f1, autoSet] = isosurface(bw_a,0.5);
    [f2, refSet]   = isosurface(bw_m,0.5);
    
    options.curvature_smoothing = 3;
    options.verb = 0;
    [Umin,Umax,Cmin,Cmax_A,Cmean_A,Cgauss_A,Normal]           = compute_curvature(autoSet,f1,options);
    
    [Umin_M,Umax_M,Cmin_M,Cmax_M,Cmean_M,Cgauss_M,Normal_M]   = compute_curvature(refSet,f2,options);
    
         clf;
     subplot(1,2,1);
     options.face_vertex_color = perform_saturation(abs(Cmean_A),1.2);
     plot_mesh(autoSet,f1, options); shading interp; colormap jet(256);
     title('Gaussian curvature A');
     subplot(1,2,2);
     %options.face_vertex_color = perform_saturation(abs(Cmin)+abs(Cmax),1.2);
     options.face_vertex_color = perform_saturation(abs(Cmean_M),1.2);
     plot_mesh(refSet, f2, options); shading interp; colormap jet(256);
     title('Gaussian curvature M');
    
    
    BA = hist(Cmax_A,-0.6:0.02:0.6);
    BB = hist(Cmax_M,-0.6:0.02:0.6);
    
    BA = BA/sum(BA);
    BB = BB/sum(BB);
    
    Bh(i) = sum(sqrt(BA.*BB))
end

save(['Bh-m'  num2str(method)],'Bh');
