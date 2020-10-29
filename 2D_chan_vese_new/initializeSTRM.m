
% Setting the figure location
set(0,'defaultfigureposition',[400 400 1000 1000]);

fPrefix = 'images';
load('DE');

for i = 1:56  % For each image set
    
    imPrefix = [fPrefix '/' num2str(i)];
    
    slices = DE{i};
    jj = 0;
    
    xS = []; yS = [];  xR = []; yR = [];
    for j = 1: size(slices,2)
        
        jj = jj +1;
        sn = slices(j);
        if(exist([imPrefix '/image' num2str(sn) 'DE.tif']))
            img = imread([imPrefix '/image' num2str(sn) 'DE.tif']);
        else
            img = imread([imPrefix '/image' num2str(sn) ' DE.tif']);            
        end
        
%         figure;imagesc(img);axis image;axis off; hold on;
        
        % First point for STRM second for region grow
        [x, y] = ginput(1);
        plot(x,y,'g*');
        
        % For STRM
        xS(jj) = x(1);
        yS(jj) = y(1);
        
        % For Region Growing
        %xR(jj) = x(1);
        %yR(jj) = y(1);
        
    end
    XS{i} = xS;
    YS{i} = yS;
    
    %XR{i} = xR;
    %YR{i} = yR;
%     close all;
end

save('XS','XS');  save('YS','YS');
%save('XR','XR');  save('YR','YR');