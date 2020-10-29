Sum = 0;

for i = 20:20
    list = dir(['./images/' num2str(i) '/' ]);
    for j = 1:size(list,1)-2
        
        flag = strcmp( list(j+2).name(end-5:end), 'DE.tif');
        if(flag ==1)
            % The counter
            Sum = Sum + 1; 
            
           
            % Reading image
            img = imread(['images/' num2str(i) '/' list(j+2).name]);
            figure;imagesc(img); axis image; axis off; colormap gray;
            indRed  = find(img(:,:,1) ==255  & img(:,:,2) ==0 & img(:,:,3) ==0);
            
            [m, n, p] = size(img);
            Contour = zeros(m, n);
            Contour(indRed) = 1;
            
            [B,L,N,A] = bwboundaries(Contour, 8, 'holes');
            [r,~] = find(A(:,N+1:end));
            [rr,~] = find(A(:,r));
            idx = setdiff(1:numel(B), [r(:);rr(:)]);
            bw = ismember(L,idx);
            
            imwrite(bw,['images/' num2str(i) '/' list(j+2).name(1:end-4)  '--bw.tif' ]);
            
        end
    end
    % DE{i}  
    close all;
end
