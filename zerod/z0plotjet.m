% visulaisation 3D des profils 0D
zs   = post.zerod;
op0d = post.z0dinput.option;
cons = post.z0dinput.cons;
geo = post.z0dinput.geo;
ts    = zs.temps;
profli = post.profil0d;
xli    = profli.xli;
t      = profli.temps;

h = findobj(0,'type','figure','tag','z0plotjet');
if isempty(h)
       h=figure('tag','z0plotjet');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
colormap('hot')


subplot(2,2,1)
switch post.z0dinput.machine
case 'JET'
  amat = interp1(zs.temps,geo.a,profli.temps,'nearest');
  rli   = profli.Raxe + amat * profli.xli;
  zplotprof(gca,t,rli,profli.nep./1e19,'color','r');
  rli   = profli.Raxe - amat * profli.xli;
  zplotprof(gca,t,rli,profli.nep./1e19,'color','r');

  liste        = {};
  liste{end+1} = 'ppf/@shot/LIDR/NE';
  liste{end+1} = 'ppf/@shot/LIDR/NEU';
  liste{end+1} = 'ppf/@shot/LIDR/NEL';   
  da           = cgcgetjet(post.z0dinput.shot,liste,'','');
  da           = da.ppf;
  zplotprof(gca,da.LIDR.NE.t,da.LIDR.NE.x,da.LIDR.NE.data./1e19, ...
            'color','b','linestyle','none','marker','o');
  zplotprof(gca,da.LIDR.NEU.t,da.LIDR.NEU.x,da.LIDR.NEU.data./1e19, ...
             'color','c','linestyle','none','marker','^');
  zplotprof(gca,da.LIDR.NEL.t,da.LIDR.NEL.x,da.LIDR.NEL.data./1e19, ...
             'color','c','linestyle','none','marker','v');

end
title(sprintf('Metis : %s@%d/density  profile', ...
             post.z0dinput.machine,post.z0dinput.shot));
ylabel('10^1^9 m^-^3')


subplot(2,2,2)
switch post.z0dinput.machine
case 'JET'
  amat = interp1(zs.temps,geo.a,profli.temps,'nearest');
  rli   = profli.Raxe + amat * profli.xli;
  zplotprof(gca,t,rli,profli.tep./1e3,'color','r');
  rli   = profli.Raxe - amat * profli.xli;
  zplotprof(gca,t,rli,profli.tep./1e39,'color','r');
  
  liste        = {};
  liste{end+1} = 'ppf/@shot/LIDR/TE';
  da          = cgcgetjet(post.z0dinput.shot,liste,'','');
  da           = da.ppf;
  zplotprof(gca,da.LIDR.TE.t,da.LIDR.TE.x,da.LIDR.TE.data./1e3, ...
            'color','m','linestyle','none','marker','o');

end
ylabel('Te (keV)');
xlabel('x (normalized radius)')

subplot(2,2,3)
switch post.z0dinput.machine
case 'JET'
  amat = interp1(zs.temps,geo.a,profli.temps,'nearest');
  rli   = profli.Raxe + amat * profli.xli;
  zplotprof(gca,t,rli,profli.tip./1e3,'color','b');	
  rli   = profli.Raxe - amat * profli.xli;
  zplotprof(gca,t,rli,profli.tip./1e3,'color','b');	
 
  liste        = {};
  liste{end+1}   = 'ppf/@shot/CXFM/TIMX';    % Ti0
  liste{end+1}   = 'ppf/@shot/CXFM/TI';      % profil de Ti
  da          = cgcgetjet(post.z0dinput.shot,liste,'','');
  da           = da.ppf;
  zplotprof(gca,da.CXFM.TI.t,da.CXFM.TI.x,da.CXFM.TI.data./1e3, ...
            'color','m','linestyle','none','marker','o');
  x0          = interp1(profli.temps,profli.Raxe(:,1),da.CXFM.TIMX.t,'nearest','extrap');
  zplotprof(gca,da.CXFM.TIMX.t,x0,da.CXFM.TIMX.data./1e3, ...
            'color','m','linestyle','none','marker','*');


end
ylabel('Ti (keV)');
xlabel('x (normalized radius)')




