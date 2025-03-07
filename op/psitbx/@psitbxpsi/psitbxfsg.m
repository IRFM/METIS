function fsg = psitbxfsg(fsd)

% PSITBXPSI/PSITBXFSG

if ~strcmp(fsd.format,'FS') | fsd.psitbxfun.grid.storage(1) ~= 'G'
 error('PSITBXFSG only applies to "FS" with "Grid" storage')
end
fact = 1;
if size(fsd.psitbxfun.grid,3) == 1, fact = 2*pi; end

% now build psitbxfun's
fsg.grho  = mean(metric(fsd,-1),[2,3],'metric');
fsg.surf  = sum(metric(fsd,[2,3]),[2,3])*fact;
fsg.area  = sum(fsd.psitbxfun.^2,2)/2;
fsg.vol   = sum(metric(fsd,'V').*fsd.psitbxfun.^2,2)*fact;
fsg.darea = diff(fsg.area,1);
fsg.dvol  = diff(fsg.vol ,1);
