function g = psitbxgrid(type,storage,x,varargin)

% PSITBXGRID   PsiTbx-Grid object constructor
% PSITBXGRID(TYPE,STORAGE,X[,PAR])

% template
if nargin == 0
 g.type = 'Undefined';
 g.label = {};
 g.storage = 'Undefined';
 g.x = {};
 [g.par,g.t] = gridpar;
 [g.csp.xk,g.csp.Mx,g.csp.Mxs]  = deal(cell(1,3));
 return
end

% type casting
if nargin == 1 & isa(type,'psitbxgrid')
 g = type;
 return
end

if nargin < 3, error('Specify TYPE, STORAGE and X'); end

% type
switch upper(type(1))
 case 'O', g.type = 'Orthogonal';        g.label = {'X','Y','Z'};
 case 'C', g.type = 'Cylindrical';       g.label = {'R','Z','Phi'};
 case 'T', g.type = 'Toroidal';          g.label = {'A','Theta','Phi'};
 case 'F', g.type = 'Flux-Coordinate';   g.label = {'Rho','Theta','Phi'};
 case 'D', g.type = 'Diagnostic-Chords'; g.label = {'s'};
 otherwise, error('Invalid coordinate TYPE')
end

% storage
switch upper(storage(1))
 case 'G', g.storage = 'Grid';
  if g.type(1) == 'D', error('"Diagnostic-Chords" cannot have "Grid" STORAGE'), end
 case 'P', g.storage = 'Points';
 case 'T', g.storage = 'Time-Points';
 otherwise, error('Invalid STORAGE')
end

% x
if ischar(x) & upper(x(1)) == 'D'%efault
 switch [g.type(1) g.storage(1)]
  case 'FG'
   x = {linspace(0,1,41),linspace(-pi,pi,129),NaN};
  case 'CG'
   x = {linspace(.614,1.147,28),linspace(-.76,.76,65)};
  otherwise
   error('No default for this TYPE and STORAGE')
 end 
end
if length(x) < 3 & g.type(1) ~= 'D', [x{end+1:3}] = deal(NaN); end
switch g.storage(1)
 case 'G'
  for k = 1:length(x), x{k} = x{k}(:)'; end
 case {'P','T'}
  dim = [];
  for k = 1:length(x)
   if length(x{k}) > 1
    if isempty(dim)
     dim = size(x{k});
    elseif ~isequal(dim,size(x{k}))
     error('For "Points" or "Time-Points" all X{:} must have the same size unless scalar')
    end
   end
  end
end
g.x = x;

% par and t
[par,t,nt] = gridpar(varargin{:});
switch g.type(1)
 case {'O','C'}
  [par.r0,par.z0,par.psi,par.fsd,par.dcd] = deal([]);
 case {'T'} % requires PAR and T
  if isempty(par.r0), error('"Toroidal" coordinates require (R0,Z0)'), end
  [par.psi,par.fsd,par.dcd] = deal([]);
 case {'F'}
  par.dcd = deal([]);
 case {'D'}
  if isempty(par.dcd), error('"Diagnostic-Chords" require a PSITBXDCD'), end
  [par.r0,par.z0,par.psi,par.fsd] = deal([]);
end
switch g.storage(1)
 case {'G','P'}, if any(g.type(1) == 'OC'), t = []; end
 case 'T'
  if nt > 1 & nt ~= dim(end)
   error('For "Time-Points" length of T must match last dimension of X')
  end
end
g.t = t; 
g.par = par;

% spline matrices
[g.csp.xk,g.csp.Mx,g.csp.Mxs]  = deal(cell(1,3));
if g.storage(1) == 'G'
 switch g.type(1)
  case 'O', ec = 'aaa';
  case 'C', ec = 'aap';
  case 'T', ec = 'spp';
  case 'F', ec = 'spp';
 end
 for k = 1:min(length(x),length(ec))
  if length(x{k}) > 1
   if ec(k) == 'p' & g.x{k}(end) - g.x{k}(1) ~= 2*pi
    ec(k) = 'a';
   end
   [g.csp.Mx{k},g.csp.xk{k}] = csdec(g.x{k},g.x{k},ec(k));
   if ec(k) == 's'
    g.csp.Mxs{k} = g.csp.Mx{k};
    g.csp.Mx{k}  = csdec(g.x{k},g.x{k},'a');
   end
  end
 end
end

g = class(g,'psitbxgrid');
