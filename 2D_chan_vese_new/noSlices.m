

Ex = [19,20,	44,	25,	27,	4,	41,	34,	1];
num = 0;
num2 = 0;
for i = 1:60
   
    if(find(Ex==i))
       num2 = num2 + size(DE{i},2);
    else
        num = num + size(DE{i},2);
    end
end
