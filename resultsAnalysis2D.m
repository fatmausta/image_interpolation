

% 3 for CMF with Bh  % 5 FWHM   % 7 Region Growing % 
% 8 for STRM+1SD
% 9 for STRM+2SD
% 10 for STRM+3SD
% 11 for STRM+4SD

method = 7;

load(['results/Dice-m' num2str(method)]);
load(['results/rmse-m' num2str(method)]);
load(['results/A_a-m'  num2str(method)]);
load('results/A_m');

for i = 1:60

    Dice(i)  = mean(sDice{i});
    rmse(i) = mean(sRmse{i});
    areaM(i)  =  mean(sA_m{i});
    areaA(i)  =  mean(sA_a{i});
    areaDiff(i) = mean(sA_a{i}-sA_m{i});
end