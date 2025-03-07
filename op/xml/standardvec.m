function [out]=standardvec(in)

%% Jo Lister September 2004
%% presents a vector or matrix in standard form [a,a,a;b,b,b]

nj = size(in,2);
ni=size(in,1);

out='[';

for i=1:ni

if(nj-1)
  out = [out sprintf('%g,',in(i,1:nj-1)) sprintf('%g;',in(i,nj))];
else 
  out=[out sprintf('%g;',in(i,1))];
end

end

out(end)=']';
