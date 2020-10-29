%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 57:61
    list = dir(['images/' num2str(i) '/']);
   jj = 0; Num = [];
   for j = 1:size(list,1)-2
        
        flag = strcmp( list(j+2).name(end-5:end), 'DE.tif');
        
        prefix = list(j+2).name(1:end-6);
        
        
        if(flag ==1)
            jj = jj + 1;
            if(prefix(1,end) == ' ')
              prefix = prefix(1,1:end-1);  
            end
            
            Num(jj) = str2double(prefix(1,6:end));
            
            
            img = imread(['images/' num2str(i) '/' prefix '.tif']);
            imgn = img(:,:,1);
            indRed  = find(img(:,:,1) ==255  & img(:,:,2) ==0 & img(:,:,3) ==0);
            indBlue = find(img(:,:,1) ==0    & img(:,:,2) ==0 & img(:,:,3) ==255);
            
            % Removing the blue colour pixels: Epicardium
            [m, n, p] = size(img);
            [y, x] = ind2sub(size(img),indBlue);
            
            epi_bw = poly2mask(x,y,m,n);
            
            for k = 1:size(indBlue,1)
                nbhood = img(y(k)-3:y(k)+3,x(k)-3:x(k)+3);
                pixelInt = median(nbhood(:));
                imgn(y(k), x(k)) = pixelInt;
            end
            
            % Removing the red colour pixels: Endocardium
            [y, x] = ind2sub(size(img),indRed);
            for k = 1:size(indRed,1)
                nbhood = img(y(k)-3:y(k)+3,x(k)-3:x(k)+3);
                pixelInt = median(nbhood(:));
                imgn(y(k), x(k)) = pixelInt;
            end
            imwrite(imgn, ['images/' num2str(i) '/' list(j+2).name(1:end-7)  '--orig.tif']);
            
        elseif(flag==0)
            % img  =  imread([num2str(i) '/Original/' list(j+2).name]);
            % figure;imagesc(img); axis image; axis off;
        end
        
    end
    DE{i} = Num;
end