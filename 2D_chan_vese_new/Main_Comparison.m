
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
        
        uu = imread(['figures/' num2str(i) '/m' num2str(method) '-' num2str(jj) '-bw.png']);
        %      h=figure;imagesc(img); hold on; axis off; axis image; colormap gray;
        %      contour(uu,[0.5 0.5],'c','LineWidth', 2);
        
        % Loading Manual segmentation
        if(exist([imPrefix '/image' num2str(sn) 'DE--bw.tif']))
            bw_m = logical(imread([imPrefix '/image' num2str(sn) 'DE--bw.tif']));
        else
            bw_m = logical(imread([imPrefix '/image' num2str(sn) ' DE--bw.tif']));
        end
        bw_a = (uu>0.5);
        %      contour(bw_m,[0.5 0.5],'y--','LineWidth', 2);
        
        % Dividing to border zone and core
        cImg = uint8(bw_m) .* img;
        Max = max(cImg(:));
        halfMax = double(Max) * 0.65;
        
        bw_m_core = cImg > halfMax;
        bw_m_gray = bw_m - bw_m_core;
        
        
        cImga = uint8(bw_a) .* img;
        Max = max(cImga(:));
        halfMax = double(Max) * 0.65;
        
        bw_a_core = cImga > halfMax;
        bw_a_gray = bw_a - bw_a_core;
        
        % Dice computation
        
        DiceTI(jj)  =  Dice_Index(bw_a,bw_m);
        DiceCI(jj)  =  Dice_Index(bw_a_core,bw_m_core);
        DiceBZ(jj)  =  Dice_Index(bw_a_gray,bw_m_gray);
        
    end
    DiceT{i} = DiceTI;
    DiceC{i} = DiceCI;
    DiceB{i} = DiceBZ;
    
    mDiceT(i) = mean(DiceTI);
    mDiceC(i) = mean(DiceCI);
    mDiceB(i) = mean(DiceBZ);
end

%Exclusions
Ex = [1,4,19,20,25,27,34,41,44];

mDiceT(Ex) = 0;
mDiceT = nonzeros(mDiceT);
mDiceC(Ex) = 0;
mDiceC = nonzeros(mDiceC);
mDiceB(Ex) = 0;
mDiceB = nonzeros(mDiceB);

display(['Mean Dice for Total Infarct: ', num2str(mean(mDiceT))]);
display(['Std Dice: ', num2str(std(mDiceT))]);

display(['Mean Dice for Core Infarct: ', num2str(mean(mDiceC))]);
display(['Std Dice: ', num2str(std(mDiceC))]);

display(['Mean Dice for Border Infarct: ', num2str(mean(mDiceB))]);
display(['Std Dice: ', num2str(std(mDiceB))]);
