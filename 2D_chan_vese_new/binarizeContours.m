%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 57:61
    
    for j = 1: 400 % Infinite loop
        
        flag = exist(['images/' num2str(i) '/image' num2str(j) '.tif']);
        
        if(flag==0)
            continue;
        end
        
        img = imread(['images/' num2str(i) '/image' num2str(j) '.tif']);
        
        imgn = img(:,:,1);
        indRed  = find(img(:,:,1) ==255  & img(:,:,2) ==0 & img(:,:,3) ==0);
        indBlue = find(img(:,:,1) ==0    & img(:,:,2) ==0 & img(:,:,3) ==255);
        
        % Removing the blue colour pixels: Epicardium
        [m, n, p] = size(img);
    
        Contour = zeros(m, n);
        Contour(indBlue) = 1;
        [B_epi, bw_epi, N, A] = bwboundaries(Contour, 8, 'holes');
    
        Contour = zeros(m, n);
        Contour(indRed) = 1;
        [B_endo, bw_endo, N, A] = bwboundaries(Contour, 8, 'holes');
    
        imwrite(bw_epi, ['images/' num2str(i) '/image' num2str(j) '--epi.tif']);
        imwrite(bw_endo,['images/' num2str(i) '/image' num2str(j) '--endo.tif']);
    end
end