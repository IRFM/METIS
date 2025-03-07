function y = metric(fsd,dim)

% PSITBXPSI/METRIC

if ~strcmp(fsd.format,'FS') | fsd.psitbxfun.grid.storage(1) ~= 'G'
 error('METRIC of a PSITBXPSI object only applies to "FS" with "Grid" storage')
end
g = fsd.psitbxfun.grid;
g = psitbxgrid(g.type,g.storage,g.x,g.t,g.par,fsd);
y = metric(g,dim);
