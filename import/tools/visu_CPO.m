function sortie=visu_CPO(shot,run,data,param,equ)
%
% sortie=visu_CPO(shot,run,[data,param])
% shot and run : UAL
% data, param, outputfile of CRONOS
%
%
%if ~exist('interpos') 
  %addpath ~basiuk/work/Projet_Cronos
  %addpath ~phuy/work/Projet_Cronos_basiuk
  %zineb_path;
%end
if nargin < 2
  run = 2;
end
if nargin < 1
  shot = 10;
end
if nargin < 5
  equ = 0;
end

idx = euitm_open ('euitm', shot, run);
disp(['reading coreprof ----------------------->'])
coreprof   = euitm_get(idx, 'coreprof');
long  = length(coreprof);
disp(['                                          length(coreprof) =',int2str(long)])
disp(['reading equilibrium ----------------------->'])
equi       = euitm_get(idx, 'equilibrium');
if isempty(equi(1).profiles_1d.gm1)
  disp(['reading equilibrium/1 ----------------------->']) 
  equi       = euitm_get(idx, 'equilibrium/1');
end
longe  = length(equi);
disp(['                                          length(equilibrium) =',int2str(longe)])
disp(['reading neoclassic ----------------------->'])
neo        = euitm_get(idx, 'neoclassic');
longb  = length(neo);
disp(['                                          length(neoclassic) =',int2str(longb)])
disp(['reading coresource ----------------------->'])
coresource = euitm_get(idx, 'coresource');
longs  = length(coresource);
disp(['                                          length(coresource) =',int2str(longs)])
disp(['reading antennas ----------------------->'])
try
   antenne    = euitm_get(idx, 'antennas');
   longa      = length(antenne);
catch 
  disp(['KO'])
  disp(['reading antennas/1 ----------------------->'])
  try 
     antenne    = euitm_get(idx, 'antennas/1');     
  catch
     disp(['KO'])
     antenne      = [];
  end
end
longa      = length(antenne);
disp(['                                          length(antennas) =',int2str(longa)])
disp(['reading coretransp ----------------------->'])
coretransp = euitm_get(idx, 'coretransp');
if isempty(coretransp(1).values(1).te_transp.diff_eff)
  disp(['reading coretransp/1 ----------------------->'])
  coretransp = euitm_get(idx, 'coretransp/1');
end
longt  = length(coretransp);
disp(['                                          length(coretransp) =',int2str(longt)])

euitm_close(idx); 
tend = coreprof(long).time;
tdeb = coreprof(1).time;
disp(['__________________________________________________________'])
disp([' ITM run from ',num2str(tdeb,4),'s to ',num2str(tend,4),' s'])
disp(['__________________________________________________________'])

clear coefee coefii
for k=1:longt
  ttr(k)      = coretransp(k).time;
  coefee(k,:) = coretransp(k).values(1).te_transp.diff_eff(:)';
  coefii(k,:) = coretransp(k).values(1).ti_transp.diff_eff(:,1)';
  coefnn(k,:) = coretransp(k).values(1).ne_transp.diff_eff(:,1)';
  coefvn(k,:) = coretransp(k).values(1).ne_transp.vconv_eff(:,1)';

end
sortie.transport.ttr = ttr;
sortie.transport.coefee = coefee;
sortie.transport.coefii = coefii;
sortie.transport.coefnn = coefnn;
sortie.transport.coefvn = coefvn;

if nargin >= 4
  indcronos    = max(iround(data.gene.temps,ttr(end)));
  indttrcronos = iround(ttr,data.gene.temps(indcronos));
  if isnan(indttrcronos)
      indttrcronos = length(ttr);
  end
else
  indttrcronos = length(ttr);
end
clear jsou xsou tsou Pesou nesou
for k=1:longs
  tsou(k)      = coresource(k).time;
end
indzer=find(diff(tsou)==0);
if ~isempty(indzer)
  tsou(indzer) = [];
  longs        = length(tsou);
end

for k=1:longs
  for l=1:length(coresource(k).values)
    jsou(k,:,l)    = coresource(k).values(l).j(:)';
    Pesou(k,:,l)   = coresource(k).values(l).qe.exp(:)';
    Pisou(k,:,l)   = coresource(k).values(l).qi.exp(:,1)';
    nesou(k,:,l)   = coresource(k).values(l).se.exp(:)';
  end
  xsou(k,:)    = coresource(k).values(1).rho_tor_norm(:)';
end
if nargin >= 4
  indcronos    = max(iround(data.gene.temps,tsou(end)));
  indsoucronos = iround(tsou,data.gene.temps(indcronos));
  if isnan(indsoucronos)
      indsoucronos = length(tsou);
  end
else
  indsoucronos = length(tsou);
end

sortie.source.tsou  = tsou;
sortie.source.jsou  = jsou;
sortie.source.Pesou = Pesou;
sortie.source.Pisou = Pisou;
sortie.source.xsou = xsou;
sortie.source.nesou  = nesou;

clear jtotsim qsim psisim times xcore Ip Itot tesim tesou tisou tisim jsim pesim
clear nesim nisim
for k=1:long
  times(k)      = coreprof(k).time;
  psisim(k,:)   = coreprof(k).psi.value(:)';
  psid1(k,:)    = coreprof(k).psi.ddrho(:)';
  psid2(k,:)    = coreprof(k).psi.d2drho2(:)';
  jsim(k,:)     = coreprof(k).psi.jni.value(:)';
  jtotsim(k,:)  = coreprof(k).profiles1d.jtot.value(:)';
  tesim(k,:)    = coreprof(k).te.value(:)';
  tebound(k,:)  = coreprof(k).te.boundary.value(1)';
  nebound(k,:)  = coreprof(k).ne.boundary.value(1)';
  tibound(k,:)  = coreprof(k).ti.boundary.value(1,1)';
  psibound(k,:) = coreprof(k).psi.boundary.value(1)';
  nesim(k,:)    = coreprof(k).ne.value(:)';
  for l = 1:size(coreprof(k).ni.value,2)
  	nisim(k,:,l)    = coreprof(k).ni.value(:,l)';
  end
  pesim(k,:)    = coreprof(k).profiles1d.pe.value(:)';
  pi1sim(k,:)    = coreprof(k).profiles1d.pi.value(:,1)';
  pitsim(k,:)    = coreprof(k).profiles1d.pi_tot.value(:)';
  densou(k,:)    = coreprof(k).ne.source_term.value(:)';
  dnedt(k,:)    = coreprof(k).ne.ddt(:)';
  dtedt(k,:)    = coreprof(k).te.ddt(:)';
  tesou(k,:)    = coreprof(k).te.source_term.value(:)';
  tisou(k,:)    = coreprof(k).ti.source_term.value(:,1)';
  tisim(k,:)    = coreprof(k).ti.value(:,1)';
  %tisim(k,:)    = coreprof(k).ti.value;
  qsim(k,:)     = coreprof(k).profiles1d.q.value(:)';
  xcore(k,:)    = coreprof(k).rho_tor_norm(:)';
  rhocore(k,:)    = coreprof(k).rho_tor(:)';
  drhodt(k,:)    = coreprof(k).drho_dt(:)';
  Itot(k)       = coreprof(k).psi.boundary.value(1);
  Ip(k)         = coreprof(k).globalparam.current_tot;
  li(k)         = coreprof(k).globalparam.li;
  betan(k)      = coreprof(k).globalparam.beta_normal;
  betap(k)      = coreprof(k).globalparam.beta_pol;
  vloop(k)      = coreprof(k).globalparam.vloop;
  amn           = coreprof(k).composition.amn;
  zn            = coreprof(k).composition.zn;
  zeff(k,:)    = coreprof(k).profiles1d.zeff.value(:)';
end
if nargin >= 4
  indcronos    = max(iround(data.gene.temps,times(end)));
  indprofcronos = iround(times,data.gene.temps(indcronos));
  if isnan(indprofcronos)
       indprofcronos = length(times);
  end
else
  indprofcronos = length(times);
end

sortie.coreprof.times   = times;
sortie.coreprof.psisim  = psisim;
sortie.coreprof.psid1sim  = psid1;
sortie.coreprof.psid2sim  = psid2;
sortie.coreprof.tisim   = tisim;
sortie.coreprof.tesim   = tesim;
sortie.coreprof.nisim   = nisim;
sortie.coreprof.nitsim   = sum(nisim,3);
sortie.coreprof.aesim    = sum(nisim,3) ./ nesim;
sortie.coreprof.nesim   = nesim;
sortie.coreprof.pesim   = pesim;
sortie.coreprof.pi1sim   = pi1sim;
sortie.coreprof.pitsim   = pitsim;
sortie.coreprof.zeff    = zeff;
sortie.coreprof.amn     = amn;
sortie.coreprof.zn      = zn;
sortie.coreprof.jsim    = jsim;
sortie.coreprof.dnedt   = dnedt;
sortie.coreprof.dtedt   = dtedt;
sortie.coreprof.jtotsim   = jtotsim;
sortie.coreprof.tesou   = tesou;
sortie.coreprof.nesou   = densou;

sortie.coreprof.tisou   = tisou;
sortie.coreprof.qsim    = qsim;
sortie.coreprof.xcore   = xcore;
sortie.coreprof.Itot    = Itot;
sortie.coreprof.Ip      = Ip;
sortie.coreprof.li      = li;
sortie.coreprof.betan   = betan;
sortie.coreprof.tebound = tebound;
sortie.coreprof.nebound = nebound;

sortie.coreprof.tibound = tibound;
sortie.coreprof.psibound = psibound;
sortie.coreprof.rho     = rhocore;
sortie.coreprof.drhodt  = drhodt;


if nargin >= 4
  indd = iround(data.gene.temps,times(1));
  indc = iround(data.gene.temps,times(end));
  tec  = data.prof.te(indc,:);
  tic  = data.prof.ti(indc,:);
  jtotc= data.prof.jmoy(indc,:);

end
clear gm1 gm2 gm3 gm4 gm5 gm6 gm7 gm8 gm9 ftr vpr qeq jeq phieq psieq rhoeq
dVdrho = zeros(longe,length(equi(1).profiles_1d.gm1));
for k=1:longe
  time(k)     = equi(k).time;
  gm1(k,:) = equi(k).profiles_1d.gm1(:)';
  gm2(k,:) = equi(k).profiles_1d.gm2(:)';
  gm3(k,:) = equi(k).profiles_1d.gm3(:)';
  gm4(k,:) = equi(k).profiles_1d.gm4(:)';
  gm5(k,:) = equi(k).profiles_1d.gm5(:)';
  gm6(k,:) = equi(k).profiles_1d.gm6(:)';
  gm7(k,:) = equi(k).profiles_1d.gm7(:)';
  gm8(k,:) = equi(k).profiles_1d.gm8(:)';
  gm9(k,:) = equi(k).profiles_1d.gm9(:)';
  ftr(k,:) = equi(k).profiles_1d.ftrap(:)'; 
  dpsidrho(k,:) = equi(k).profiles_1d.dpsidrho_tor(:)';
  vpr(k,:) = equi(k).profiles_1d.vprime(:)';   
  spr(k,:) = equi(k).profiles_1d.aprime(:)';   
  qeq(k,:) = equi(k).profiles_1d.q(:)';
  jeq(k,:) = equi(k).profiles_1d.jphi(:)';
  phieq(k,:) = equi(k).profiles_1d.phi(:)';
  psieq(k,:) = equi(k).profiles_1d.psi(:)';
  rhoeq(k,:) = equi(k).profiles_1d.rho_tor(:)';
    volume(k,:) = equi(k).profiles_1d.volume(:)';
  if isempty( equi(k).profiles_1d.dvdrho)
   dVdr = z0polyderive(cat(2,-1,rhoeq(k,:)),cat(2,0,volume(k,:)),7);
   %%%[yout,dVdr] = interpos(rhoeq(k,:),volume(k,:),-0.1,[1 2],[0 volume(k,end)]);
   dVdrho(k,:) = dVdr;
  else
    dVdrho(k,:) =  equi(k).profiles_1d.dvdrho(:)';
  end
  Fdia(k,:)  = equi(k).profiles_1d.F_dia(:)';
end

sortie.equi.times      = time;
sortie.equi.gm1        = gm1;
sortie.equi.gm2        = gm2;
sortie.equi.gm3        = gm3;
sortie.equi.gm4        = gm4;
sortie.equi.gm5        = gm5;
sortie.equi.gm6        = gm6;
sortie.equi.gm7        = gm7;
sortie.equi.gm8        = gm8;
sortie.equi.gm9        = gm9;
sortie.equi.vpr        = vpr;
sortie.equi.spr        = spr;
sortie.equi.jeq        = jeq;
sortie.equi.qeq        = qeq;
sortie.equi.dpsidrho   = dpsidrho;
sortie.equi.phi        = phieq;
sortie.equi.rho        = rhoeq;
sortie.equi.psi        = psieq;
sortie.equi.volume     = volume;
sortie.equi.dVdrho     = dVdrho;
sortie.equi.F          = Fdia;
grid1=equi(longe).profiles_2d(1).grid.dim1;
grid2=equi(longe).profiles_2d(1).grid.dim2;
psieq = equi(longe).profiles_2d(1).psi;

clear jboot sig timeb neoee
for k=1:longb
  timeb(k)     = neo(k).time;
  jboot(k,:)   = neo(k).jboot(:)';
  sig(k,:)     = neo(k).sigma(:)';
  xneo(k,:)    = xcore(k,:);
  neoee(k,:)   = neo(k).te_neo.diff_eff(:)';
  if ~isempty(neo(k).ti_neo.diff_eff)
    neoii(k,:)   = neo(k).ti_neo.diff_eff(:,1)';
    neonn(k,:)   = neo(k).ne_neo.diff_eff(:)';
    neovn(k,:)   = neo(k).ne_neo.vconv_eff(:)';
  else
    neoii(k,:) =    neoee(k,:) * 0;
    neonn(k,:) =    neoee(k,:) * 0;
    neovn(k,:) =    neoee(k,:) * 0;
  end
end
sortie.neo.times = timeb;
sortie.neo.jboot = jboot;
sortie.neo.sigma = sig;
sortie.neo.xneo  = xneo;
sortie.neo.coefee = neoee;
sortie.neo.coefnn = neonn;
sortie.neo.coefvn = neovn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = findobj(0,'type','figure','tag','cpo');
if isempty(h)
       h=figure('tag','cpo');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

subplot(3,3,1)

if nargin >= 4
  plot(param.gene.x(1:5:end),data.prof.psi(indc,1:5:end)*2*pi,'r',xcore(indprofcronos,:),psisim(indprofcronos,:),'bo')
  legend(['CRONOS '],['ETS '])
  title(['psi at t=',num2str(times(indprofcronos),4),' s'])
  A       =  data.prof.psi(indc,:)*2*pi;
  B       =  psisim(indprofcronos,:);
  [errvp, ierrvp] = max(abs(A-B)./abs(A));
  text('units','normalized','position',[0.1 0.9],'string',[' error = ',num2str(100*errvp,3),' %, index=',int2str(ierrvp)]);

else
  plot(xcore(indprofcronos,:),psisim(indprofcronos,:),'k')
  title(['psi at t=',num2str(times(indprofcronos),4),' s'])

end
ylabel('A/m²')

subplot(3,3,2)

if nargin >= 4
  plot(times,psisim(:,10),'r',data.gene.temps(indd:indc),data.prof.psi(indd:indc,10)*2*pi,'bo')
  legend('ETS','CRONOS')
  title('Evolution of psi at x=0.1')
else
  plot(times,psisim(:,10),'r')  
  title('Evolution of psi at x=0.1')
end

ylabel('A/m²')
xlabel('t (s)')

subplot(3,3,3)
if nargin >= 4
  plot(times,psisim(:,90),'r',data.gene.temps(indd:indc),data.prof.psi(indd:indc,90)*2*pi,'bo')
  legend('ETS','CRONOS')
  title('Evolution of psi at x=0.9')
else
  ind90 = min(find(rhoeq(1,:) > 0.9));
  plot(times,psisim(:,ind90),'r')  
  title('Evolution of psi at x=0.9')
end

ylabel('A/m²')
xlabel('t (s)')

subplot(3,3,4)

newspr=spr.*dpsidrho.*rhocore(1,end);
for k=1:longb

  Iboot(k)=trapz(xneo(k,:),jboot(k,:).*newspr(end,:))*rhocore(end,end);

end
% -> to be validated______________
%for k=1:long
%  Ip(k)  = Ip(k) * rhocore(k,end)
%end
%_________________________________
if nargin >= 4
  if long > 2
    plot(times(2:end),li(2:end),'r',times(2:end),betap(2:end),'b',times(2:end),Ip(2:end)/1e6,'g',timeb,Iboot/1e6,'b')
  else
    plot(times(2:end),li(2:end),'rs',times(2:end),betap(2:end),'bs',times(2:end),Ip(2:end)/1e6,'gs',timeb,Iboot/1e6,'bs')
  end
  hold on
  temps=data.gene.temps(indd:indc);
  lic = data.gene.li(indd:indc)/data.gene.li(4)*data.equi.li(4);
  betapc = data.gene.betap(indd:indc);
  Ipc= data.cons.ip(indd:indc)/1e6;
  Ibootc=data.gene.iboot(indd:indc)/1e6;
  plot(temps,lic,'ro',temps,betapc,'bo',temps,Ipc,'g+',temps,Ibootc,'k*')
  hold off
  legend('li','b_p','I_p','I_b_o_o_t')
  title('Evolution of some average quantities (line ETS, mark CRONOS)')
else
plot(times(2:end),li(2:end),times(2:end),betap(2:end),times(2:end),Ip(2:end)/1e6,timeb,Iboot/1e6)
legend('li','b_p','I_p','I_b_o_o_t')
end
xlabel('t (s)')

subplot(3,3,5)

if nargin >= 4
  plot(xcore(indprofcronos,:),jtotsim(indprofcronos,:),'r',param.gene.x(1:5:end),jtotc(1:5:end),'bo')
  title(['Profile of <J> at t=',num2str(times(indprofcronos),4),' s'])
  legend(['CRONOS '],['ETS '])
else
  plot(xcore(1,:),jtotsim(1,:),'ro',xcore(end,:),jtotsim(end,:),'r')
  legend(['t=',num2str(times(1),4),' s'],['t=',num2str(times(end),4),' s'])
  title('Profile of <J>')
end
ylabel('A/m²')
xlabel('x')

subplot(3,3,6)

if nargin >= 4
  plot(times(:),jtotsim(:,10),'r',data.gene.temps(indd:indc),data.prof.jmoy(indd:indc,10),'bo')  
  legend(['ETS '],['CRONOS '])
else
  plot(times(:),jtotsim(:,10),'r')  
end
title('Evolution of <J> at x=0.1')
ylabel('A/m²')
xlabel('t (s)')

subplot(3,3,7)



if nargin >= 4 
  plot(xcore(indprofcronos,:),qsim(indprofcronos,:),'r',param.gene.x(1:5:end),data.prof.q(indc,1:5:end),'bo')
  legend('ETS','CRONOS')
  title(['q at t=',num2str(times(indprofcronos),4),' s'])
else
  plot(xcore(1,:),qsim(1,:),'ro',xcore(end,:),qsim(end,:),'k')
  legend(['t=',num2str(times(1),4),' s'],['t=',num2str(times(end),4),' s'])
  title('q')
end
xlabel('x')

subplot(3,3,8)

newvpr = vpr .* dpsidrho;
if nargin >= 4
  plot(param.gene.x,data.equi.vpr(indc,:),'r',xcore(end,1:10:end),newvpr(end,1:10:end),'bo',xcore(end,:),dVdrho(end,:),'b--')
  title(['dV/drho at t=',num2str(times(end),4),' s'])
  legend('CRONOS','ETS','ETS, from volume')
  A       =  data.equi.vpr(indc,2:end);
  B       =   newvpr(end,2:end);
  [errvp, ierrvp] = max(abs(A-B)./A);
  text('units','normalized','position',[0.1 0.9],'string',[' error = ',num2str(100*errvp,3),' %, index=',int2str(ierrvp)]);

else
  plot(xcore(1,:),newvpr(1,:),'r',xcore(end,:),newvpr(end,:),'bo',xcore(end,:),dVdrho(end,:),'b--')
  legend(['t=',num2str(times(1),4),' s'],['t=',num2str(times(end),4),' s'],['t=',num2str(times(end),4),' s'])
  title('dV/drho')
end

subplot(3,3,9)

if ~isempty(psieq)
if (~(any(find(size(psieq) == 1))))
    if equ == 0
       contourf(grid1,grid2,psieq')
    else
       contourf(grid1,grid2,psieq)
    end
end
end
ylabel('Z(m)');
xlabel('R(m)');
axis('square')
if nargin >= 4
  hold on
  if nargin == 5
    plot(data.geo.R(indc,:),data.geo.Z(indc,:),'w-')
  else
    plot(data.geo.R(indc,:),data.geo.Z(indc,:)-data.geo.z0(indc),'w-')
  end
  title('last equilibrium, white: CRONOS')
else
  title('last equilibrium')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = findobj(0,'type','figure','tag','kinetic');
if isempty(h)
       h=figure('tag','kinetic');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

subplot(4,2,1)

plot(xcore(indprofcronos,:),tesim(indprofcronos,:),'r')

if nargin >= 4

  hold on
  plot(param.gene.x,data.prof.te(indc,:),'bo')
  legend('ETS','CRONOS')
  A  =  tesim(indprofcronos,:);
  B  =  data.prof.te(indc,:);
  [errvp, ierrvp] = max(abs(A-B)./A);
  text('units','normalized','position',[0.1 0.9],'string',[' error = ',num2str(100*errvp,3),' %, index=',int2str(ierrvp)]);

end
title(['Te profile at t=',num2str(times(indprofcronos),4),' s'])

ylabel('eV')

subplot(4,2,2)

plot(xcore(indprofcronos,:),tisim(indprofcronos,:),'r')

if nargin >= 4

  hold on
  plot(param.gene.x,data.prof.ti(indc,:),'bo')
  legend('ETS','CRONOS')
  A  =  tisim(indprofcronos,:);
  B  =  data.prof.ti(indc,:);
  [errvp, ierrvp] = max(abs(A-B)./A);
  text('units','normalized','position',[0.1 0.9],'string',[' error = ',num2str(100*errvp,3),' %, index=',int2str(ierrvp)]);
end
title(['Ti profile at t=',num2str(times(indprofcronos),4),' s'])

ylabel('eV')

subplot(4,2,3)
ind10 = min(find(rhoeq(1,:)>0.1));
ind50 = min(find(rhoeq(1,:)>0.5));

if nargin >= 4
  plot(times,tesim(:,ind10),'r',data.gene.temps(indd:indc),data.prof.te(indd:indc,10),'bo')
  legend('ETS','CRONOS')
  title('Evolution of Te at x=0.1')
  %A  =  tesim(:,10);
  %B  =  data.prof.te(indd:indc,10);
  %[errvp, ierrvp] = max(abs(A-B)./A);
  %text('units','normalized','position',[0.1 0.9],'string',[' error = ',num2str(100*errvp,3),' %, index=',int2str(ierrvp)]);
else
  plot(times,tesim(:,ind10),times,tesim(:,ind50))
  legend('x=0.1','x=0.5')
  title('Evolution of Te')
end
ylabel('eV')
xlabel('t (s)')

subplot(4,2,4)

if nargin >= 4
  plot(times,tisim(:,10),'r',data.gene.temps(indd:indc),data.prof.ti(indd:indc,10),'bo')
  legend('ETS','CRONOS')
  title('Evolution of Ti at x=0.1')
  %A  =  tisim(:,10);
  %B  =  data.prof.ti(indd:indc,10);
  %[errvp, ierrvp] = max(abs(A-B)./A);
  %text('units','normalized','position',[0.1 0.9],'string',[' error = ',num2str(100*errvp,3),' %, index=',int2str(ierrvp)]);
else
  plot(times,tisim(:,1),times,tisim(:,ind50))
  legend('x=0.1','x=0.5')
  title('Evolution of Ti')
end

ylabel('eV')
xlabel('t (s)')

subplot(4,2,5)

plot(xcore(indprofcronos,:),nesim(indprofcronos,:),'r')

if nargin >= 4

  hold on
  plot(param.gene.x,data.prof.ne(indc,:),'bo')
  legend('ETS','CRONOS')
  A  =  nesim(indprofcronos,:);
  B  =  data.prof.ne(indc,:);
  [errvp, ierrvp] = max(abs(A-B)./A);
  text('units','normalized','position',[0.1 0.9],'string',[' error = ',num2str(100*errvp,3),' %, index=',int2str(ierrvp)]);
end

title(['ne profile at t=',num2str(times(indprofcronos),4),' s'])

ylabel('m^-^3')

subplot(4,2,6)

plot(xcore(indprofcronos,:),sum(nisim(indprofcronos,:,:),3),'r')

if nargin >= 4

  hold on
  plot(param.gene.x,data.prof.ni(indc,:),'bo')
  legend('ETS','CRONOS')
  A  =  sum(nisim(indprofcronos,:,:),3);
  B  =  data.prof.ni(indc,:);
  [errvp, ierrvp] = max(abs(A-B)./A);
  text('units','normalized','position',[0.1 0.9],'string',[' error = ',num2str(100*errvp,3),' %, index=',int2str(ierrvp)]);
end
title(['sum(ni) profile at t=',num2str(times(indprofcronos),4),' s'])

ylabel('m^-^3')


subplot(4,2,7)

if nargin >= 4
  plot(times,nesim(:,10),'r',data.gene.temps(indd:indc),data.prof.ne(indd:indc,10),'bo')
  legend('ETS','CRONOS')
  title('Evolution of ne at x=0.1')
else
  plot(times,nesim(:,ind10),times,nesim(:,ind50))
  legend('x=0.10','x=0.5')
  title('Te')
end
ylabel('m^-^3')
xlabel('t (s)')


subplot(4,2,8)

if nargin >= 4
  plot(times,sum(nisim(:,10,:),3),'r',data.gene.temps(indd:indc),data.prof.ni(indd:indc,10,:),'bo')
  legend('ETS','CRONOS')
  title('Evolution of sum(ni) at x=0.1')
else
  plot(times,sum(nisim(:,ind10,:),3))
  legend('x=0.1','x=0.5')
  title('Te')
end
ylabel('m^-^3')
xlabel('t (s)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = findobj(0,'type','figure','tag','cpo1');
if isempty(h)
       h=figure('tag','cpo1');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

plot(xneo(end,:),jboot(end,:)/1e6,'r--',xcore(end,:),jtotsim(end,:)/1e6,'k')
if nargin >= 4

  hold on
  plot(param.gene.x, data.neo.jboot(indc,:)/1e6,'ro',...
       param.gene.x, data.prof.jmoy(indc,:)/1e6,'ko')
  hold off
  legend('ETS J_b_o_o_t','ETS (jtot/value), <J>','CRONOS J_b_o_o_t','CRONOS <J>')
  B = jtotsim(end,1:end-5)/1e6;
  A = data.prof.jmoy(indc,1:end-5)/1e6;
  [errvp, ierrvp] = max(abs(A-B)./A);
  text('units','normalized','position',[0.1 0.9],'string',[' error_J = ',num2str(100*errvp,3),' %, index=',int2str(ierrvp)]);


else
  legend('ETS J_b_o_o_t','ETS (jtot/value), <J>')
end
ylabel('MA/m^2')
title(['current profiles at t=',num2str(times(end),4),' s'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = findobj(0,'type','figure','tag','cpo2');
if isempty(h)
       h=figure('tag','cpo2');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
plot(times(1:end-1),diff(times)*1e3)
if nargin >= 4
  hold on
  plot(data.gene.temps(indd:indc),data.gene.dt(indd:indc)*1e3,'ro')
  legend('ETS','CRONOS')
  hold off
end
title('dt evolution')
xlabel('s')
ylabel ('ms')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = findobj(0,'type','figure','tag','coef');
if isempty(h)
       h=figure('tag','coef');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

subplot(1,2,1)
if nargin >= 4
  plot(xcore(indttrcronos,:),coefee(indttrcronos,:),'r',param.gene.x,data.coef.ee(indc,:)./data.prof.ne(indc,:),'ro', ...
       xcore(indttrcronos,:),coefii(indttrcronos,:),'b',param.gene.x,data.coef.ii(indc,:)./data.prof.ni(indc,:),'bo',...
       xcore(indttrcronos,:),coefnn(indttrcronos,:),'m',param.gene.x,data.coef.nn(indc,:),'mo')
   legend('ee_I_T_M','ee_C_R_O','ii_I_T_M','ii_C_R_O','nn_I_T_M','nn_C_R_O')
else
   plot(xcore(end,:),coefee(end,:),xcore(end,:),coefii(end,:),xcore(end,:),coefnn(end,:))
   legend('ee_I_T_M','ii_I_T_M','nn_I_T_M')
end
title(['Anormal transport coefficient, t=',num2str(ttr(indttrcronos),4),' s'])
ylabel('m²/s')

subplot(1,2,2)
if nargin >= 4
  plot(xcore(indttrcronos,:),coefvn(indttrcronos,:),'r',param.gene.x,data.coef.vn(indc,:),'ro')
   legend('vn_I_T_M','vn_C_R_O_N_O_S')
else
   plot(xcore(end,:),coefvn(end,:))
   legend('vn_I_T_M')
end
title(['pinch, t=',num2str(ttr(indttrcronos),4),' s'])
ylabel('m/s')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = findobj(0,'type','figure','tag','source');
if isempty(h)
       h=figure('tag','source');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

plot(xcore(indsoucronos,:),tesou(indsoucronos,:),'b',xcore(indsoucronos,:),tisou(indsoucronos,:),'r')
if nargin >= 4
  hold on
  plot(param.gene.x,data.source.totale.el(indc,:),'bo',param.gene.x,data.source.totale.ion(indc,:),'ro')

  legend('electrons (ETS)','ions (ETS)','electrons (CRONOS)','ions (CRONOS)')
  title(['Total heating source profiles at t=',num2str(tsou(indsoucronos),4),' s'])
  B = tesou(indsoucronos,:);
  A = data.source.totale.el(indc,:);
  maxa=max(abs(A));
  [errvp, ierrvp] = max(abs(A-B)/maxa);
  text('units','normalized','position',[0.2 0.9],'string',[' error_e_l_e_c_t_r_o_n = ',num2str(100*errvp,3),' %, index=',int2str(ierrvp)]);
  B = tisou(indsoucronos,:);
  A = data.source.totale.ion(indc,:);
  maxa=max(abs(A));
  [errvp, ierrvp] = max(abs(A-B)/maxa);
  text('units','normalized','position',[0.2 0.85],'string',[' error_i_o_n = ',num2str(100*errvp,3),' %, index=',int2str(ierrvp)]);
else
  legend('te','ti')
  title('heat source')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = findobj(0,'type','figure','tag','source_par');
if isempty(h)
       h=figure('tag','source_par');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
if nargin >= 4
    indcronos=max(find(data.gene.temps<tsou(end)));
    inditm=iround(tsou,data.gene.temps(indcronos));
    if isnan(inditm)
        inditm = length(tsou);
    end
else
   inditm = length(tsou); 
end
    
subplot(3,3,1)

plot(xsou(inditm,:),Pesou(inditm,:,1),'r')
title('Heat sources, ECRH')

if nargin >= 4
 indg = 1:5:length(param.gene.x);
  hold on
  plot(param.gene.x(indg),data.source.fce.el(indcronos,indg),'ro')
  title('Heat sources, ECRH, line:ETS, mark:CRONOS, blue:electron, red:ion')
  hold off
end

subplot(3,3,2)

plot(xsou(inditm,:),Pesou(inditm,:,2),'r',xsou(inditm,:),Pisou(inditm,:,2),'b')
title('ICRH')

if nargin >= 4
  hold on
  plot(param.gene.x(indg),data.source.fci.el(indcronos,indg),'ro',param.gene.x(indg),data.source.fci.ion(indcronos,indg),'bo')
  hold off
end

subplot(3,3,3)

plot(xsou(inditm,:),Pesou(inditm,:,3),'r')
title('LH')

if nargin >= 4
  hold on
  plot(param.gene.x(indg),data.source.hyb.el(indcronos,indg),'ro')
  hold off
end

subplot(3,3,4)

plot(xsou(inditm,:),Pesou(inditm,:,4),'r',xsou(inditm,:),Pisou(inditm,:,4),'b')
title('NBI')

if nargin >= 4
  hold on
  plot(param.gene.x(indg),data.source.idn.el(indcronos,indg),'ro',param.gene.x(indg),data.source.idn.ion(indcronos,indg),'bo')
  hold off
end

subplot(3,3,5)

plot(xsou(inditm,:),Pesou(inditm,:,6),'r',xsou(end,:),Pisou(inditm,:,6),'b')
title('Cold neutral')

if nargin >= 4
  hold on
  plot(param.gene.x(indg),data.source.n0.el(indcronos,indg),'ro',param.gene.x(indg),data.source.n0.ion(indcronos,indg),'bo')
  hold off
end

subplot(3,3,6)

plot(xsou(inditm,:),Pesou(inditm,:,7),'r')
title('equipartition')
if nargin >= 4
  hold on
  plot(param.gene.x(indg),data.source.qei(indcronos,indg),'ro')
  hold off
end

subplot(3,3,7)

plot(xsou(inditm,:),Pesou(inditm,:,8),'r')
title('ohmic power')
if nargin >= 4
  hold on
  plot(param.gene.x(indg),data.source.ohm(indcronos,indg),'ro')
  hold off
end

subplot(3,3,8)

plot(xsou(inditm,:),Pesou(inditm,:,9),'r')
title('Bremsstrahlung')
if nargin >= 4
  hold on
  plot(param.gene.x(indg),data.source.brem(indcronos,indg),'ro')
  hold off
end

subplot(3,3,9)

plot(xsou(inditm,:),Pesou(inditm,:,10),'r')
title('line radiation')
if nargin >= 4
  hold on
  plot(param.gene.x(indg),data.source.prad(indcronos,indg),'ro')
  hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = findobj(0,'type','figure','tag','source_par2');
if isempty(h)
       h=figure('tag','source_par2');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])


subplot(2,2,1)

plot(xsou(inditm,:),nesou(inditm,:,6),'r')
title('Cold neutral density source')

if nargin >= 4
  hold on
  plot(param.gene.x(indg),data.source.n0.ne(indcronos,indg),'ro')
  hold off
end


subplot(2,2,2)

plot(xsou(inditm,:),nesou(inditm,:,4),'r')
title('NBI density source')

if nargin >= 4
  hold on
  plot(param.gene.x(indg),data.source.idn.ne(indcronos,indg),'ro')
  hold off
end

subplot(2,2,3)

plot(xcore(inditm,:),densou(inditm,:),'r')
title('totale density source')

if nargin >= 4
  hold on
  plot(param.gene.x(indg),data.source.totale.ne(indcronos,indg),'ro')
  hold off
end


subplot(2,2,4)

plot(times,nebound,'r',times,nesim(:,end),'k--')
title('electron density')
legend('ETS boundary','ETS edge value')
if nargin >= 4
  hold on
  plot(data.gene.temps(indd:indc),data.cons.ne1(indd:indc),'ro',data.gene.temps(indd:indc),data.prof.ne(indd:indc,end),'k*')
  hold off
  legend('ETS boundary','ETS edge value','CRONOS boundary','CRONOS edge value')
end


nefromion=zeros(size(xcore));
for l=1:length(times)
  for k=1:length(zn)
    nefromion(l,:) = nefromion(l,:) + zn(k) * nisim(l,:,k);
  end  
end
niqei = zeros(size(xcore));
for l=1:length(times)
  for k=1:length(zn)
    niqei(l,:) = niqei(l,:) + zn(k).^2 * nisim(l,:,k) ./ amn(k);
  end  
end
lnei  = 14.9 - 0.5 * log(nesim/1e20) + log(tesim/1e3);
epsi0 =  8.854187817620391e-12;
e     =  1.602176462000000e-19;
me    = 9.109381880000001e-31;
mp    = 1.672648500000000e-27;
taue  = (12 * pi.^1.5 / sqrt(2)) .* (epsi0.^2 ./ e.^4 .* sqrt(me)) .* ((e.*tesim).^1.5 ./ nesim ./ lnei);
qei   = 3 .* me ./ mp ./ taue .* niqei .* (e .* tesim-e .* tisim(:,:,1));
if nargin >= 4

  for l=1:length(data.gene.temps)
    niqeicro(l,1:length(param.gene.x)) = 0.;
    for k=1:length(zn)
      niqeicro(l,:) = niqeicro(l,:) + param.compo.z(k).^2 * data.impur.impur(l,:,k) ./ param.compo.a(k);
    end  
  end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if any(abs(nefromion ./nesim - 1) > 1e-10)
  h = findobj(0,'type','figure','tag','neutrality');
if isempty(h)
       h=figure('tag','neutrality');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
  
  plot(nefromion ./nesim)
  title('Electroneutrality for the whole profile')
  
else
  disp('neutrality OK')  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = findobj(0,'type','figure','tag','boundary');
if isempty(h)
       h=figure('tag','boundary');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

subplot(2,2,1)

plot(times,rhocore(:,end),'ro')

if nargin >= 4

  hold on
  plot(data.gene.temps(indd:indc),data.equi.rhomax(indd:indc),'k')
  legend('ITM','CRONOS')

else
  legend('ITM')
end
title('rhomax')
ylabel('m')

subplot(2,2,2)

plot(times,psisim(:,end),'ro')

if nargin >= 4

  hold on
  plot(data.gene.temps(indd:indc),data.prof.psi(indd:indc,end)*2*pi,'k')
  legend('ITM','CRONOS')

else
  legend('ITM')
end
hold off
title('psi edge')

subplot(2,2,3)

plot(times,tebound,'ro')

if nargin >= 4

  hold on
  if any(data.mode.pe(indd:indc) == 2)
    plot(data.gene.temps(indd:indc),data.prof.te(indd:indc,end),'k')
  
      
  else
    plot(data.gene.temps(indd:indc),data.cons.te1(indd:indc),'k')
  
  end
  legend('ITM','CRONOS')

else
  legend('ITM')
end
hold off
title('Te boundary')
ylabel('eV')

subplot(2,2,4)
pi2         = pi*pi;
const       = - 8 * pi2 * 1e-7;

psid1_1_ITM = const * psibound ./ (dVdrho(:,end).*gm2(:,end));
plot(times,psid1_1_ITM,'ro')

if nargin >= 4
  psid1_1_CRO = const * data.cons.ip ./ (data.equi.c2c(:,end));

  hold on
  plot(data.gene.temps(indd:indc),psid1_1_CRO(indd:indc,end),'k')
  legend('ITM','CRONOS')

else
  legend('ITM')
end
hold off
title('psi boundary (for Ip)')





