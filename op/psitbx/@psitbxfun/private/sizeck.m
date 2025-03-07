function [dimf,ks] = sizeck(dimf,dimg)

ks = [];
for k = 1:length(dimg)
 if dimg(k) == 1 & (k > length(dimf) | dimf(k) ~= dimg(k))
  ks = [ks k];
  dimf = [dimf(1:k-1) 1 dimf(k:end)];
 elseif dimf(k) ~= dimg(k)
  dimf = [];
  return
 end
end
