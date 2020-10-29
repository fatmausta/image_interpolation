

num = 0;
for i= 2:22
   TS(i) =  sum(Time{i});
   num = num + size(Time{i},2);
end