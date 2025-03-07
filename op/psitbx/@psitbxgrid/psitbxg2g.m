function g2 = psitbxg2g(g1,g2type,varargin)

% PSITBXG2G   PsiTbx-Grid transformation
% PSITBXG2G(GRID,TYPE)

% psitbxg2g(gC,'T',{r0,z0}|gT)
% psitbxg2g(gT,'C'[,{r0,z0}|gT])
% psitbxg2g(gC,'F',{r0,z0,fC}|psi01)

switch upper([g1.type(1),g2type(1)])

 case {'OO','CC','TT','FF','DD'}, g2 = g1;

 case {'OD','CD','TD','FD'}
  error('Transformation to "Diagnostic-Chords" not possible')

 case {'OT','OF','TO','TF','FO','DO','DT','DF'}
  g2 = psitbxg2g(psitbxg2g(g1,'C',varargin{:}),g2type,varargin{:});
 case 'FC'
  g2 = psitbxg2g(psitbxg2g(g1,'T',varargin{:}),g2type,varargin{:});

 case 'OC'
  [par,t,nt] = gridpar(varargin{:},g1);
  [x,g2storage] = g2p(g1.x,g1.storage);
  g2 = psitbxgrid(g2type,g2storage,...
   {sqrt(x{1}.^2+x{2}.^2),x{3},atan2(x{2},x{1}),x{4:end}},t);

 case 'CO'
  [par,t,nt] = gridpar(varargin{:},g1);
  [x,g2storage] = g2p(g1.x,g1.storage);
  if all(isnan(x{3}))
   x = {x{1},repmat(NaN,size(x{2})),x{2},x{4:end}};
  else
   x = {x{1}.*cos(x{3}),x{1}.*sin(x{3}),x{2},x{4:end}};
  end
  g2 = psitbxgrid(g2type,g2storage,x,t);

 case 'CT'
  [par,t,nt] = gridpar(varargin{:},g1);
  [x,g2storage] = g2t(g1.x,g1.storage,nt);
  k0 = ndims(x{1});
  g2 = psitbxgrid(g2type,g2storage,...
   {sqrt(xop(par.r0',k0,'-',x{1},':').^2+xop(par.z0',k0,'-',x{2},':').^2),...
    atan2(xop(x{2},':','-',par.z0',k0),xop(par.r0',k0,'-',x{1},':')),...
    x{3:end}},par,t);

 case 'TC'
  [par,t,nt] = gridpar(varargin{:},g1);
  [x,g2storage] = g2t(g1.x,g1.storage,nt);
  k0 = ndims(x{1});
  g2 = psitbxgrid(g2type,g2storage,...
   {xop(par.r0',k0,'-',x{1}.*cos(x{2}),':'),...
    xop(par.z0',k0,'+',x{1}.*sin(x{2}),':'),...
    x{3:end}},t);

 case 'CF'
  [par,t,nt] = gridpar(varargin{:},g1);
  [x,g2storage] = g2t(g1.x,g1.storage,nt);
  k0 = ndims(x{1});
  psi1 = psitbxf2f(par.psi,g1);
  g2 = psitbxgrid(g2type,g2storage,...
   {sqrt(reshape(psi1.x,size(x{1}))),...
    atan2(xop(x{2},':','-',par.psi.zmag',k0),xop(par.psi.rmag',k0,'-',x{1},':')),...
    x{3:end}},par,t);
    
 case 'FT'
  [par,t,nt] = gridpar(varargin{:},g1);
  [x,g2storage] = g2t(g1.x,g1.storage,nt);
  fsd1 = psitbxf2f(par.fsd,g1);
  g2 = psitbxgrid(g2type,g2storage,...
   {reshape(fsd1.x,size(x{1})),...
    x{2:end}},par,t);

 case 'DC'
  [par,t,nt] = gridpar(varargin{:},g1);
  dcd = par.dcd; t = dcd.t;
  x = g1.x; g2storage = g1.storage;
		nd = size(x{1},1);
  x = {rep(dcd.rd,nd) - x{1}.*rep(cos(dcd.pvd).*cos(dcd.tvd),nd),...
       rep(dcd.zd,nd) + x{1}.*rep(sin(dcd.pvd),nd),...
       x{1}.*rep(cos(dcd.pvd).*sin(dcd.tvd),nd)};
  x = {sqrt(x{1}.^2+x{3}.^2),x{2},mod(atan2(x{3},x{1})+rep(dcd.phid,nd),2*pi)};
  g2 = psitbxgrid(g2type,g2storage,x,t);  

 otherwise, error('Invalid new grid type')
end

function x = rep(x,n)
x = repmat(reshape(x,[1,size(x)]),[n,1]);
