function rmse = RMSE(manualCont,autoCont)
%This function computes the Mean Absolute distance between a manual and
%a automatic contour
%

for k= 1:size(manualCont,1)
    d(:,k) = sqrt((autoCont(:,1)-manualCont(k,1)).^2+(autoCont(:,2)-manualCont(k,2)).^2);
end 
minDist = min(d);
sqError = minDist.^2;
MeansqError = mean(sqError);
rmse = sqrt(MeansqError);

end