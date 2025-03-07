function f = scd(op,f,dim,dometric,pad)

% to perform sum, cumsum and diff

if nargin < 4, dometric = 'no'; end
if nargin < 5, pad = 0; end
dometric = strcmp(lower(dometric),'metric');

if f.grid.storage(1) ~= 'G' & f.grid.type(1) ~= 'D'
 error('DIFF, SUM and CUMSUM only possible for "Grid" storage or "Diagnostic-Chords"')
end
x = f.grid.x;
dimx = size(f.grid);

if dometric, f = f .* metric(f.grid,dim); end
f = reshape(f);
dimz = size(f); ndimz = length(dimz);

if isequal(dim,'t'), dim = ndimz; end
if any(dim > ndimz), error('Invalid dimension'), end
dim(dimz(dim) == 1) = []; % not worth working on singleton
dim = fliplr(sort(dim));
kper = 1:ndimz; kper(dim) = [];
z = permute(f.x,[dim,kper]);
if op == 's', z(isnan(z)) = 0; end

for k = 1:length(dim)
 dimz = size(z);
 d = [];
 if op == 's' | op == 'c'
  if dim(k) == ndimz & length(f.t) > 1 % along time
   d = (f.t([2:end,end])-f.t([1,1:end-1]))'/2;
   if op == 's', f.t = []; end
  elseif dim(k) > length(dimx) % along a non grid
  elseif (f.grid.type(1) == 'D')
   d = (x{dim(k)}([2:end,end],1) - x{dim(k)}([1,1:end-1],1)) / 2;
  else
   d = (x{dim(k)}([2:end,end])   - x{dim(k)}([1,1:end-1]))'  / 2;
  end
 end
 if op == 's'
  if dim(k) <= 3
   x{dim(k)} = NaN;
  else
   x(dim(k)) = [];
  end
 end
 z = feval(['scd' op],z,d,pad);
 z = permute(z,[2:ndimz 1]);
end

if op == 's'
 f.x = squeeze(z);
 f.grid = psitbxgrid(f.grid.type,f.grid.storage,x,f.grid.par,f.grid.t);
else
 f.x = squeeze(ipermute(z,[kper,dim]));
end

function z = scds(z,d,pad) % sum
if ~isempty(d)
 dimz = size(z);
 z = z .* repmat(d,[1,dimz(2:end)]);
end
z = sum(z,1);
function z = scdc(z,d,pad) % cumsum
if ~isempty(d)
 dimz = size(z);
 z = z .* repmat(d,[1,dimz(2:end)]);
end
z = cumsum(z,1);
function z = scdd(z,d,pad) % diff
dimz = size(z);
z = reshape([repmat(pad,1,prod(dimz(2:end)));z(2:end,:)-z(1:end-1,:)],dimz);
