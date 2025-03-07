function psi = psitbxpsi(x,g,t,form,rzpsimag,iter,tol)

% PSITBXPSI   PsiTbx Poloidal Flux object constructor
% PSITBXPSI(X,GRID,T,FORM[,[RMAG;ZMAG;PSIMAG],ITER,TOL])

% type casting
if nargin == 1 & isa(x,'psitbxpsi')
 psi = x;
 return
end

if nargin < 1, x = []; end
if nargin < 2, g = psitbxgrid; end
if nargin < 3, t = []; end
if nargin < 4, form = 'Undefined'; end
if nargin < 5, rzpsimag = []; end
if nargin < 6, iter = 20; end
if nargin < 7, tol = 1e-4; end

if size(rzpsimag,1) >= 1, rmag = rzpsimag(1,:); else, rmag = []; end
if size(rzpsimag,1) >= 2, zmag = rzpsimag(2,:); else, zmag = []; end
if size(rzpsimag,1) >= 3, pmag = rzpsimag(3,:); else, pmag = []; end

switch form
 case {'10','+0','-0','01','FS','Undefined'}
 otherwise, error('Invalid psitbxpsi object format')
end

psi.format = form;
psi.rmag   = rmag;
psi.zmag   = zmag;
psi.psimag = pmag;
psi.iter   = iter;
psi.tol    = tol;
psi = class(psi,'psitbxpsi',psitbxfun(x,g,t));
