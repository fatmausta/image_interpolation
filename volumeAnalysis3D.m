addpath('imMinkowski');

% 3 for CMF with Bh
% 5 FWHM
% 7 Region Growing %
% 8 for STRM+1SD
% 9 for STRM+2SD
% 10 for STRM+3SD
% 11 for STRM+4SD

M = [3,5 ,7,8,9,10,11];

for j = 1:7
    
    method = M(j);
    for i =   1:60
        
        bw_a =     metaImageRead(['results/'  num2str(i) '/m'  num2str(method)  '-bw_a.mhd']);
        bw_m =     metaImageRead(['results/'  num2str(i) '/'  '-bw_m.mhd']);
        
        bw_a = logical(bw_a);
        bw_m = logical(bw_m);
        
        volume_a(i) = nnz(bw_a(:));
        volume_m(i) = nnz(bw_m(:));
    end
    
    
    save(['volume-m' ],'volume_m');
    save(['volume-m' num2str(method)],'volume_a');
end