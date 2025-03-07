function psitbxput(ver,shot,redo)

% PSITBXPUT		Computes and stores data in the RESULTS tree
%   PSITBXPUT(VER,SHOT[,REDO]) User should not invoke this function. It is
%   called automatically when accessing an empty or not up to date node
%   \RESULTS::PSITBX:*

if nargin < 3, redo = 1; end

if 1
 if isempty(mdsopen('RESULTS',shot)), return, end
 if ~redo & ver == mdsdata('\RESULTS::PSITBX:AS:VERSION_NUM'), return, end
 psi = psitbxtcv(shot,(0:100)/20,'+0');
 mdssetdef('\RESULTS::PSITBX')
else
 if isempty(mdsopen('pcs',shot)), return, end
 psi = psitbxtcv(shot,'RAMP');
 mdssetdef('\PCS::MGAMS.PSITBX')
end
if isempty(psi), return, end
fsd = psitbxp2p(psi,'FS');
rzs = psitbxg2g(fsd.grid,'C',fsd);
fsg = psitbxfsg(fsd);

put('rmag'  ,fsd.rmag	,{psi.t}			    ,'m')
put('zmag'  ,fsd.zmag	,{psi.t}			    ,'m')
put('psimag',fsd.psimag ,{psi.t}			    ,'Vs')
put('as'    ,fsd.x	,{fsd.grid.x{1},fsd.grid.x{2},psi.t},'m')
put('rs'    ,rzs.x{1}	,{fsd.grid.x{1},fsd.grid.x{2},psi.t},'m')
put('zs'    ,rzs.x{2}	,{fsd.grid.x{1},fsd.grid.x{2},psi.t},'m')
put('grho'  ,fsg.grho.x ,{fsd.grid.x{1} 	     ,psi.t},'m^-1')
put('surf'  ,fsg.surf.x ,{fsd.grid.x{1} 	     ,psi.t},'m^2')
put('area'  ,fsg.area.x ,{fsd.grid.x{1} 	     ,psi.t},'m^2')
put('darea' ,fsg.darea.x,{fsd.grid.x{1} 	     ,psi.t},'m^2')
put('vol'   ,fsg.vol.x  ,{fsd.grid.x{1} 	     ,psi.t},'m^3')
put('dvol'  ,fsg.dvol.x ,{fsd.grid.x{1} 	     ,psi.t},'m^3')

mdsclose

function put(node,x,dim,unit)
n = length(dim);
expr = ['BUILD_SIGNAL(BUILD_WITH_UNITS($1,$2),'];
if n > 1, expr = [expr sprintf(',$%d',3:n+1)]; end
expr = [expr sprintf(',BUILD_WITH_UNITS($%d,"s"))',n+2)];
if 1
 mdsput([node ':FOO'],expr,'x',mdscvt(x,'f'),unit,dim{:})
 mdsput([node ':VERSION_NUM'],psitbxver,'f')
else
 mdsput(node,expr,'x',mdscvt(x,'f'),unit,dim{:})
 mdsput('VER',ver,'f')
end
