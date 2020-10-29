fPrefix = 'images';
load('DE');
options.GPU = 0;
set(0,'defaultfigureposition',[400 400 600 600]);

method = 7;
% Method Select
% 1 Chan-Vese method
% 2 Log Likelihood method
% 3 CMF with Bh
% 4 Half-Max Threshold Myo
% 5 FWHM   
% 6 Non-Seg Myo+1SD Threshold Seg
% 7 Region Growing
% 8 STRM+1SD
% 9 STRM+2SD
% 10 STRM+3SD
% 11 STRM+4SD

Filt = [0,3,5,10,15,30,40,50,100,200,300];

for k = 1:size(Filt,2)
    
    for i = 1:60
        
        % For each image set
        imPrefix = [fPrefix '/' num2str(i)];
        
        slices = DE{i};
        jj = 0;
        
        Dice = [];
        for j = 1: size(slices,2)
            
            jj = jj +1;
            sn = slices(j);
            
            img = imread([imPrefix '/image' num2str(sn) '--orig.tif']);
            IMG = imread([imPrefix '/image' num2str(sn) '.tif']);
            
            % Loading Manual segmentation
            if(exist([imPrefix '/image' num2str(sn) 'DE--bw.tif']))
                bw_m = logical(imread([imPrefix '/image' num2str(sn) 'DE--bw.tif']));
            else
                bw_m = logical(imread([imPrefix '/image' num2str(sn) ' DE--bw.tif']));
            end
            
            uu = imread(['figures/' num2str(i) '/m' num2str(method) '-' num2str(jj) '-bw.png']);
            
%             h = figure; imagesc(img); hold on; axis off; axis image; colormap gray;
%             contour(uu,[0.5 0.5],'c','LineWidth', 2);
%             contour(bw_m,[0.5 0.5],'y--','LineWidth', 2);

            bw_a = (uu>0.5);  
            
            % Connectivity Filtering
            bw_a = bwareaopen(bw_a, Filt(k));
            
            % Dice computation
            DiceTI(jj)  =  Dice_Index(bw_a,bw_m);
            
        end
        
        DiceT{i} = DiceTI;        
        mDiceT(i) = mean(DiceTI);
        
    end
    
    %Exclusions
    Ex = [1,4,19,20,25,27,34,41,44];
    
    mDiceT(Ex) = 0;
    mDiceT = nonzeros(mDiceT);
    Needed(k) = mean(mDiceT);
    
    display(['Mean Dice for Total Infarct (Filt=', num2str(Filt(k)), '): ', num2str(mean(mDiceT))]);
    display(['Std Dice: ', num2str(std(mDiceT))]);
end
