function f = metric(g,dim)

% PSITBXGRID/METRIC

if g.type(1) == 'D', f = 1; return, end
if g.storage(1) ~= 'G', error('Only for "Grid" storage'), end
if strcmp(dim,'V'), cphi = [1/2,1/3]; dim = 3; else, cphi = [1,1]; end
dim = sort(dim); dim(dim > 3) = [];
f = 1;
[xx,tmp] = g2t(g.x(1:3),g.storage,length(g.t));
kt = ndims(xx{1});
switch g.type(1)
 case 'O'
 case 'C'
  if any(dim == 3)
   f = xx{1};
  elseif dim == -3
   f = 1./xx{1};
  end
 case 'T'
  if any(dim == 2)
   f = xx{1};
  end
  if any(dim == 3)
   f = f .* xop(cphi(1)*g.par.r0',kt,'-',cphi(2)*xx{1}.*cos(xx{2}),':');
  end
  if dim == -2
   f = 1./xx{1};
  elseif dim == -3
   f = 1 ./ xop(cphi(1)*g.par.r0',kt,'-',cphi(2)*xx{1}.*cos(xx{2}),':');
  end
 case 'F'
  if isempty(g.par.fsd) | isempty(g.par.r0)
   error('Requires "fsd" and "r0" parameters for "Flux" coordinates')
  end
  a  = g.par.fsd.x;
  ar = getfield(psitbxf2f(g.par.fsd,[1,0,0]),'x');
  at = getfield(psitbxf2f(g.par.fsd,[0,1,0]),'x');
  at = sqrt(a.^2 + at.^2);
  r  = xop(cphi(1)*g.par.r0',kt,'-',cphi(2)*a.*cos(xx{2}),':');
  warning off MATLAB:divideByZero
  switch sprintf('%1d',dim)
   case '1'  , f = ar;
   case '2'  , f = at;
   case '3'  , f = r;
   case '12' , f = a .* ar;
   case '13' , f = ar .* r;
   case '23' , f = at .* r;
   case '123', f = a .* ar .* r;
   case '-1' , f = at ./ a ./ ar; k = find(a == 0); f(k) = 1./ar(k);
   case '-2' , f = 1 ./ a;
   case '-3' , f = 1 ./ r;
  end
  warning on MATLAB:divideByZero
end
if ~isequal(f,1)
 f = psitbxfun(f,g,g.t);
end
