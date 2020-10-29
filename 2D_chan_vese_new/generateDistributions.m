fPrefix = 'images';

% Image stat

load('DE');

for i = 2:40  % For each image set
    
    imPrefix = [fPrefix '/' num2str(i)];
    
    slices = DE{i};
    jj = 0;
    intScar = [];
    for j = 1: size(slices,2)
        
        jj = jj +1;
        sn = slices(j);
        
        img = imread([imPrefix '/image' num2str(sn) '--orig.tif']);
        
        % Scar region
        if(exist([imPrefix '/image' num2str(sn) 'DE--bw.tif']))
            bw_scar = imread([imPrefix '/image' num2str(sn) 'DE--bw.tif']);
        else
            bw_scar = imread([imPrefix '/image' num2str(sn) ' DE--bw.tif']);
        end
            
        indx = (bw_scar ==1);
        intScar{jj} = img(indx);
        
        % Normal myocardium region
        bw_epi  = imread([imPrefix '/image' num2str(sn) '--epi.tif']);
        bw_endo = imread([imPrefix '/image' num2str(sn) '--endo.tif']);
        
        bw_myo = bw_epi - bw_endo;
        bw_myo = logical(bw_myo);
        bw_normal = bw_myo - bw_scar;
        indx = (bw_normal ==1);
        intNormal{jj} = img(indx);
        
    end
    IntScar{i}   = cell2mat(intScar');
    IntNormal{i} = cell2mat(intNormal');
end

IntsScar = cell2mat(IntScar(2:end)');
IntsNormal = cell2mat(IntNormal(2:end)');

samples = 0:1:255;

intPDF_scar   = ksdensity(double(IntsScar),samples,'Width',3);
intPDF_normal = ksdensity(double(IntsNormal),samples,'Width',3);

IntsTotal = [IntsScar; IntsNormal];

intPDF_total = ksdensity(double(IntsTotal),samples,'Width',3);

save('intPDF_scar','intPDF_scar');
save('intPDF_normal','intPDF_normal');
save('intPDF_total','intPDF_total');

figure; plot(intPDF_scar); hold on; plot(intPDF_normal,'r');
plot(intPDF_total,'k'); hold off;

