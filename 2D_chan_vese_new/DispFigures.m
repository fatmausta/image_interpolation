%

M = [3,5,7,8,9,10];

for k = 1:6
    method = M(k);
    for i = 1:5
        img = imread(['26/m' num2str(method) '-' num2str(i) '-img.png']);
        img = img(345:577,495:755,1:3);
        imwrite(img,['figures/26_m' num2str(method) '-' num2str(i) '-img.png']);
    end
end