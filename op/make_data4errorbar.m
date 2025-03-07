function y = make_data4errorbar(xm,xp)

ss = size(xm);
y = NaN * ones(ss(1),3 * ss(2));
if nargin > 1
  vn = NaN * ones(ss(1),ss(2));
else
  vn = xm;
  xp = xm;
end
for k = 1:ss(1)
  u        = cat(1,xm(k,:),xp(k,:),vn(k,:));
  y(k,:)   = u(:)';
end 