
fPrefix = 'images';
load('DE');
options.GPU = 0;

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

M = [3,5,7,8,9,10];

for k = 1:6
    disp(['Method: ', num2str(M(k))])    
    method = M(k);

for i = 26:26  % For each image set
    
    imPrefix = [fPrefix '/' num2str(i)];
    
    slices = DE{i};
    jj = 0;
    
    Dice = [];  A_m = [];  A_a = [];  rmse = [];
    for j = 1: size(slices,2)
        
        jj = jj +1;
        sn = slices(j);
        
        img = imread([imPrefix '/image' num2str(sn) '--orig.tif']);
        
        bw_epi  = imread([imPrefix '/image' num2str(sn) '--epi.tif']);
        bw_endo = imread([imPrefix '/image' num2str(sn) '--endo.tif']);
        
        bw_myo = bw_epi - bw_endo(:,:,1);
        
        indx = img(logical(bw_myo));
        Std = std(double(indx));
        Mean = mean(double(indx));
        
        Thres = Mean;
        
        mask = bw_myo;
        mask(img< Thres) =0;
        
        if(method == 1)
            % Chan-Vese method
            [uu, uu_sq, mv,stime, erriter, num] = SegmentationCV(double(img), bw_myo, mask);
            
        elseif(method == 2)
            % Log Likelihood method
            load('intPDF_scar');
            load('intPDF_normal');
            [uu, uu_sq,stime, erriter, num] = SegmentationLog2(double(img), bw_myo, mask, intPDF_scar, intPDF_normal);
            
        elseif(method == 3)
            %
            load('intPDF_scar');
            load('intPDF_normal');
            [uu, uu_sq,stime, erriter, num] = SegmentationBh(double(img), bw_myo, mask, intPDF_scar, intPDF_normal,options);
            
        elseif(method == 4)
            img_bw = double(img) .* double(bw_myo)/255;
            thres1 = max(img_bw(:))/2;
            uu = img_bw > thres1;
            
        elseif(method == 5)
            if(exist([imPrefix '/image' num2str(sn) 'DE--bw.tif']))
                bw_m = logical(imread([imPrefix '/image' num2str(sn) 'DE--bw.tif']));
            else
                bw_m = logical(imread([imPrefix '/image' num2str(sn) ' DE--bw.tif']));
            end
            img_bw = double(img) .* double(bw_m);
            thres1 = max(img_bw(:))/2;
            
            Img_bw = double(img) .* double(bw_myo)/255;
            uu = Img_bw > thres1;
            
        elseif(method == 6)
            if(exist([imPrefix '/image' num2str(sn) 'DE--bw.tif']))
                bw_m = logical(imread([imPrefix '/image' num2str(sn) 'DE--bw.tif']));
            else
                bw_m = logical(imread([imPrefix '/image' num2str(sn) ' DE--bw.tif']));
            end
            img_bw = double(img) .* double(bw_m);
            region = double(bw_myo)/255 - double(bw_m);
            img_reg = region .* double(img);
            
            Int = median(nonzeros(img_reg(:)));
            Std = std(nonzeros(img_reg(:)));
            thres2 = Int + 1 *Std;
            uu = img_bw > thres2;
            
        elseif(method == 7)
            if(exist([imPrefix '/image' num2str(sn) 'DE--bw.tif']))
                bw_m = logical(imread([imPrefix '/image' num2str(sn) 'DE--bw.tif']));
            else
                bw_m = logical(imread([imPrefix '/image' num2str(sn) ' DE--bw.tif']));
            end
            
            % Loading the points for Region Growing
            load('XR');   load('YR');
            x = XR{i}(jj);
            y = YR{i}(jj);
            
            bw_myo = double(bw_myo)/255;
            img = double(img)/255;
            imgS =  (1- (bw_myo)) * 1e4;
            imgS((bw_myo==1)) = img((bw_myo==1)) ;
            
            % Region Growing
            % Finding the seed point for region growing
            uu = region_growing(imgS,round(y),round(x),0.1);
            % uu = regiongrowing(double(img)/255,[round(y) round(x)]);
            
        elseif(method == 8)
            if(exist([imPrefix '/image' num2str(sn) 'DE--bw.tif']))
                bw_m = logical(imread([imPrefix '/image' num2str(sn) 'DE--bw.tif']));
            else
                bw_m = logical(imread([imPrefix '/image' num2str(sn) ' DE--bw.tif']));
            end
            
            img_bw = double(img) .* double(bw_myo)/255;
            
            % Loading the points for STRM
            load('XS');   load('YS');
            x = round(XS{i}(jj));
            y = round(YS{i}(jj));
            
            roi = double(img(y-4:y+4,x-4:x+4));
            Im   = mean(roi(:));
            Istd = std(roi(:));
            
            Thres = Im + (Istd * 1);
            uu = img_bw > Thres;
            
        elseif(method == 9)
            if(exist([imPrefix '/image' num2str(sn) 'DE--bw.tif']))
                bw_m = logical(imread([imPrefix '/image' num2str(sn) 'DE--bw.tif']));
            else
                bw_m = logical(imread([imPrefix '/image' num2str(sn) ' DE--bw.tif']));
            end
            
            img_bw = double(img) .* double(bw_myo)/255;
            
            
            % Loading the points for STRM
            load('XS');   load('YS');
            x = round(XS{i}(jj));
            y = round(YS{i}(jj));
            
            roi = double(img(y-4:y+4,x-4:x+4));
            Im   = mean(roi(:));
            Istd = std(roi(:));
            
            Thres = Im + (Istd * 4);
            uu = img_bw > Thres;
            
        elseif(method == 10)
            if(exist([imPrefix '/image' num2str(sn) 'DE--bw.tif']))
                bw_m = logical(imread([imPrefix '/image' num2str(sn) 'DE--bw.tif']));
            else
                bw_m = logical(imread([imPrefix '/image' num2str(sn) ' DE--bw.tif']));
            end
            
            img_bw = double(img) .* double(bw_myo)/255;
            
            
            % Loading the points for STRM
            load('XS');   load('YS');
            x = round(XS{i}(jj));
            y = round(YS{i}(jj));
            
            roi = double(img(y-4:y+4,x-4:x+4));
            Im   = mean(roi(:));
            Istd = std(roi(:));
            
            Thres = Im + (Istd * 5);
            uu = img_bw > Thres;
            
        elseif(method == 11)
            if(exist([imPrefix '/image' num2str(sn) 'DE--bw.tif']))
                bw_m = logical(imread([imPrefix '/image' num2str(sn) 'DE--bw.tif']));
            else
                bw_m = logical(imread([imPrefix '/image' num2str(sn) ' DE--bw.tif']));
            end
            
            img_bw = double(img) .* double(bw_myo)/255;
            
            
            % Loading the points for STRM
            load('XS');   load('YS');
            x = XS{i}(jj);
            y = YS{i}(jj);
            
            roi = double(img(y-4:y+4,x-4:x+4));
            Im   = mean(roi(:));
            Istd = std(roi(:));
            
            Thres = Im + (Istd * 4);
            uu = img_bw > Thres;
        end
        
%         h = figure; imagesc(img); hold on; axis off; axis image; colormap gray;
%         contour(uu,[0.5 0.5],'c','LineWidth', 2);
        
        % Loading Manual segmentation
        if(exist([imPrefix '/image' num2str(sn) 'DE--bw.tif']))
            bw_m = logical(imread([imPrefix '/image' num2str(sn) 'DE--bw.tif']));
        else
            bw_m = logical(imread([imPrefix '/image' num2str(sn) ' DE--bw.tif']));
        end
        bw_a = (uu>0.5);

%         contour(bw_m,'y:','LineWidth', 1); hold off;
 
        % Metrics%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Dice Coefficient
        Dice(jj) = Dice_Index(bw_m,bw_a);
        
        % Area computation
        A_m(jj) = nnz(bw_m);
        A_a(jj) = nnz(bw_a);
        
        % Contours of the segmentation
        C_m   = bwboundaries(bw_m);
        C_a   = bwboundaries(bw_a);
        
        % Root Mean Square Error
        if(isempty(C_a) == 1)
            rmse(jj) = 10;
        else
            rmse(jj) = RMSE(C_a{1},C_m{1});
        end
        % Saving figures
%         saveas(h, ['results/' num2str(i) '/m' num2str(method) '-' num2str(jj)  '-img.png']);
%         imwrite(bw_a,['results/' num2str(i) '/m' num2str(method) '-' num2str(jj) '-bw.png' ]);
%         imwrite(bw_m,['results/' num2str(i) '/man'  '-' num2str(jj) '-bw.png'              ]);
    end
    % Computing the mean Dice
    mDice(i) = mean(Dice);
    
    % Metric for each subject
    sDice{i}      =   Dice;
    sRmse{i}    =   rmse;
    sA_m{i}      =   A_m;
    sA_a{i}       =   A_a;
    
    %  pause(3);
    close all;
end

% Saving the  data
% save(['results/' 'Dice-m' num2str(method)],'sDice');
% save(['results/' 'rmse-m' num2str(method)],'sRmse');
% save(['results/' 'A_a-m'  num2str(method)],'sA_a');
% save(['results/' 'A_m'],'sA_m');

display(['Mean Dice: ', num2str(mean(nonzeros(mDice)))]);
display(['Std Dice: ', num2str(std(nonzeros(mDice)))]);

end