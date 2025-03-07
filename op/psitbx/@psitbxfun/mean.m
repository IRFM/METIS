function f = mean(f,dim,w,dometric)

% PSITBXFUN/MEAN

domet = 0;
switch nargin
 case 2
  w = 1;
 case 3
  if ischar(w)
   domet = strcmp(upper(w),'METRIC');
   w = 1;
  end
 case 4
  domet = strcmp(upper(dometric),'METRIC')
end
 
if ~isa(w,'double'), w = psitbxf2f(w,f.grid); end
if domet, w = metric(f.grid,dim) .* w; end
f = sum(f.*w,dim) ./ sum(~isnan(f).*w,dim);
