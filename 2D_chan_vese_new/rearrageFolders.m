% Rearrange folders

for i = 24:56
    movefile(['images/' num2str(i) '/Original/*'], ['images/' num2str(i) '/']);
end