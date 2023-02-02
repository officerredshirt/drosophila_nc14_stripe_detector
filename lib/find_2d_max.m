function [val,ix] = find_2d_max(img)

[rowmax,rowix] = max(img,[],1);
[val,colix] = max(rowmax);

ix = [rowix(colix),colix];

end