function [par,t,nt] = gridpar(varargin)

[par.r0,par.z0,par.psi,par.fsd,par.dcd,t,nt] = deal([]);
for k = nargin:-1:1
 v = varargin{k};
 switch class(v)
  case 'psitbxgrid'
   t = ifempty(v.t,t);
   v = v.par;
  case 'psitbxfun'
   t = ifempty(v.t,t);
   v = v.grid.par;
  case 'psitbxpsi'
   t = ifempty(v.t,t);
   par.r0 = ifempty(v.rmag,par.r0);
   par.z0 = ifempty(v.zmag,par.z0);
   if strcmp(v.format,'01'), par.psi = v; end
   if strcmp(v.format,'FS'), par.fsd = v; end
   v = v.grid.par;
  case 'psitbxdcd'
   t = ifempty(v.t,t);
   par.dcd = v;
 end
 switch class(v)
  case 'double'
   t = ifempty(v,t);
  case 'cell'
   if length(v) == 2
    par.r0 = ifempty(v{1},par.r0);
    par.z0 = ifempty(v{2},par.z0);
   end
  case 'struct'
   if isfield(v,'r0') , par.r0  = ifempty(v.r0 ,par.r0 ); end
   if isfield(v,'z0') , par.z0  = ifempty(v.z0 ,par.z0 ); end
   if isfield(v,'psi'), par.psi = ifempty(v.psi,par.psi); end
   if isfield(v,'fsd'), par.fsd = ifempty(v.fsd,par.fsd); end
   if isfield(v,'dcd'), par.dcd = ifempty(v.dcd,par.dcd); end
 end
end
t = t(:)'; nt = length(t);
par.r0 = par.r0(:)';
par.z0 = par.z0(:)';

if length(par.r0) ~= length(par.z0) | (length(par.r0) > 1 & length(par.r0) ~= nt)
 error('Invalid (R0,Z0) found in parameters')
end

%%%
function x = ifempty(x1,x)
if ~isempty(x1), x = x1; end
