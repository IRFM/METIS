% script du  plot des densites 0d
% compatibilite ascendante
if ~isfield(post.zerod,'nbar_nat') || all(~isfinite(post.zerod.nbar_nat))
     ulh = 0.25;
     fh  = post.zerod.modeh;
     post.zerod.nbar_nat   = min(post.zerod.negr,max(1e13, 1e20 .* (post.z0dinput.geo.b0 ./ post.zerod.q95 ./ post.z0dinput.geo.R) .^ 0.6  .*  (ulh + (1 - ulh) .* fh)));
     post.zerod.frac_pellet = post.zerod.frac_pellet .* (post.zerod.frac_pellet > 1e-2);
end


zs   = post.zerod;
cons = post.z0dinput.cons;
geo  = post.z0dinput.geo;
exp0d  = post.z0dinput.exp0d;
t    = zs.temps;
nhesup = zs.esup_fus ./ (3.56e6 .*1.6022e-19).*2 ./ zs.vp;


% precalcul
if post.z0dinput.mode_exp == 0
   vpr  = data.equi.vpr;
   nvpr = trapz(param.gene.x,vpr,2);

   ind1 = find(param.compo.z ==1);
   if isempty(ind1)
      n1m = zeros(size(t));
   else
      n1m = trapz(param.gene.x,vpr .* squeeze(sum(data.impur.impur(:,:,ind1),3)),2)./ nvpr;
   end

   indD = find(param.compo.a ==2);
   if isempty(indD)
      nDm = zeros(size(t));
   else
      nDm = trapz(param.gene.x,vpr .* squeeze(sum(data.impur.impur(:,:,indD),3)),2)./ nvpr;
   end

   indT = find(param.compo.a == 3 & param.compo.z == 1);
   if isempty(indT)
      nTm = zeros(size(t));
   else
      nTm = trapz(param.gene.x,vpr .* squeeze(sum(data.impur.impur(:,:,indT),3)),2)./ nvpr;
   end

   indhe = find(param.compo.z == 2);
   if isempty(indhe)
      nhem = zeros(size(t));
   else
      nhem = trapz(param.gene.x,vpr .* squeeze(sum(data.impur.impur(:,:,indhe),3)),2)./ nvpr;
   end

   indimp = find(param.compo.z > 2);
   if isempty(indimp)
      nimpm = zeros(size(t));
   else
      nimpm = trapz(param.gene.x,vpr .* squeeze(sum(data.impur.impur(:,:,indimp),3)),2)./ nvpr;
   end
end


% low density limit due to runaway electron instability
% ref: C. Paz-Soldan et al, NF 56 (2016) 056010
phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)
ntrap = trapz(post.profil0d.xli,post.profil0d.nep .* post.profil0d.ftrap .* post.profil0d.vpr,2) ./  ...
	max(1,trapz(post.profil0d.xli,post.profil0d.nep .* post.profil0d.vpr,2));
% expression Rosenbluth : correction du champ critique
corp = 1 + interp1(post.profil0d.temps,ntrap,zs.temps,'linear','extrap') ./ zs.nem;
lng         = 14.9 - 0.5 .* log(zs.nem ./ 1e20) + log(zs.tem ./ 1e3);
epar        = max(eps,abs(zs.vloop) ./ 2 ./ pi ./ geo.R);
vth         = min(phys.c,sqrt(zs.tem .* phys.e ./ phys.me));
ec0         = corp  .* phys.e .^ 3 .* lng ./ 4 ./ pi ./ phys.epsi0 .^ 2 ./ phys.me ./ vth .^ 2;
nlow_crit   = min(max(zs.negr),zs.nbar ./ zs.nem .* epar ./ ec0 ./ 0.07);

h = findobj(0,'type','figure','tag','z0plotn');
if isempty(h)
       h=figure('tag','z0plotn');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

k    = 9;
if post.z0dinput.mode_exp == 0

   subplot(k,1,1);	
   plot(t,zs.nbar/1e19,'r',data.gene.temps,data.gene.nbar ./ 1e19,'-.b',t,zs.negr ./ 1e19,'g',t,zs.nbar_nat ./ 1e19,'g-.',t,zs.nsat ./ 1e19,'c',t,real(cons.nbar)./1e19,'m:',t,nlow_crit ./ 1e19,'k-.');
   ylabel('nbar (1e19 m^-^3)')
   %legend('n_{bar}','n_{bar cronos}','n_{Greenwald}','n_{natural}', 'n_{OH saturation}', 'n_{reference}','Location','best','orientation','horizontal')
   %legend(gca,'boxoff')

   title(sprintf('Zerod : %s@%d/density ( r -> 0D, b: -> cronos)', ...
             post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,-t,'r',t,-t,'-.b',t,-t,'g',t,-t,'g-.',t,-t,'c',t,-t,'c:',t,-t,'k-.');
   legend('n_{bar}','n_{bar experiment}','n_{Greenwald}','n_{natural}', 'n_{OH saturation}', 'n_{reference}','n_{low limit RE instability}','Location','south','orientation','horizontal')
   hold on
   plot(t,zs.nem./1e19,'r',data.gene.temps,data.gene.nemoy./1e19,':b');
   set(gca,'ylim',[0,Inf]);
   ylabel('<Ne> (1e19 m^-^3)');
   subplot(k,1,3);	
   plot(t,zs.ne0./1e19,'r',data.gene.temps,data.prof.ne(:,1)./1e19,':b');
   ylabel('Ne0 (1e19 m^-^3)');
   subplot(k,1,4);	
   plot(t,zs.nebord./1e19,'r',data.gene.temps,data.prof.ne(:,end)./1e19,':b');
   ylabel('Nea (1e19 m^-^3)');
   subplot(k,1,5);	
   plot(t,zs.n1m./1e19,'r',data.gene.temps,n1m./1e19,':b');
   ylabel('<N_H_D_T> (1e19 m^-^3)');
   subplot(k,1,6);	
   plot(t,zs.nDm./1e19,'r',data.gene.temps,nDm./1e19,':b');
   ylabel('<N_D> (1e19 m^-^3)');
   subplot(k,1,7);	
   plot(t,zs.nTm./1e19,'r',data.gene.temps,nTm./1e19,':b');
   ylabel('<N_T> (1e19 m^-^3)');
   subplot(k,1,8);	
   plot(t,zs.nhem./1e19,'r',data.gene.temps,nhem./1e19,':b',t,nhesup./1e19,'m');
   ylabel('<N_H_e> (1e19 m^-^3)');
   subplot(k,1,9);	
   plot(t,zs.nimpm./1e19,'r',data.gene.temps,nimpm./1e19,':b');
   ylabel('<Nimp> (1e19 m^-^3)');
   xlabel('time (s)');
else

   subplot(k,1,1);	
   plot(t,zs.nbar/1e19,'r',t,exp0d.nbar ./ 1e19,'-.b',t,zs.negr ./ 1e19,'g',t,zs.nbar_nat ./ 1e19,'g-.',t,zs.nsat ./ 1e19,'c',t,real(cons.nbar)./1e19,'c:',t,nlow_crit ./ 1e19,'k-.');
   ylabel('nbar (1e19 m^-^3)')
   %legend('n_{bar}','n_{bar experiment}','n_{Greenwald}','n_{natural}', 'n_{OH saturation}', 'n_{reference}','Location','south','orientation','horizontal')
   %legend(gca,'boxoff')
   title(sprintf('Zerod : %s@%d/density ( r -> 0D, b: -> experiment)', ...
             post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,-t,'r',t,-t,'-.b',t,-t,'g',t,-t,'g-.',t,-t,'c',t,-t,'c:',t,-t,'k-.');
   legend('n_{bar}','n_{bar experiment}','n_{Greenwald}','n_{natural}', 'n_{OH saturation}', 'n_{reference}','n_{low limit RE instability}','Location','south','orientation','horizontal')
   hold on
   plot(t,zs.nem./1e19,'r',t,exp0d.nem./1e19,':b');
   set(gca,'ylim',[0,Inf]);
   ylabel('<Ne> (1e19 m^-^3)');
   subplot(k,1,3);	
   plot(t,zs.ne0./1e19,'r',t,exp0d.ne0./1e19,':b');
   ylabel('Ne0 (1e19 m^-^3)');
   subplot(k,1,4);	
   plot(t,zs.nebord./1e19,'r',t,exp0d.nebord./1e19,':b');
   ylabel('Nea (1e19 m^-^3)');
   subplot(k,1,5);	
   plot(t,zs.n1m./1e19,'r',t,exp0d.n1m./1e19,':b');
   ylabel('N_H_D_T (1e19 m^-^3)');
   subplot(k,1,6);	
   plot(t,zs.nDm./1e19,'r',t,exp0d.nDm./1e19,':b');
   ylabel('N_D (1e19 m^-^3)');
   subplot(k,1,7);	
   plot(t,zs.nTm./1e19,'r',t,exp0d.nTm./1e19,':b');
   ylabel('N_T (1e19 m^-^3)');
   subplot(k,1,8);	
   %plot(t,zs.nhem./1e19,'r',t,exp0d.nhem./1e19,':b',t,nhesup./1e19,'m');
   if isfield(zs, 'nhe4m') && ~all(zs.nhe4m(:) == 0)
           plot(t,zs.nhem./1e19,'r',t,exp0d.nhem./1e19,':b',t,nhesup./1e19,'m',t,zs.nhe4m./1e19,'g');
           legend('n_{he3}','n_{he experiment}','n_{hesup}','n_{he4m}', 'Location','northwest','orientation','horizontal')
   else
           plot(t,zs.nhem./1e19,'r',t,exp0d.nhem./1e19,':b',t,nhesup./1e19,'m');
   end
   ylabel('N_H_e (1e19 m^-^3)');
   subplot(k,1,9);	
   plot(t,zs.nimpm./1e19,'r',t,exp0d.nimpm./1e19,':b');
   ylabel('Nimp (1e19 m^-^3)');
   xlabel('time (s)');

end
hh = findobj(h,'type','axes');
for l=1:k
   lim = get(hh(l),'ylim');
   set(hh(l),'ylim',[0,max(max(lim),0.1)]);
end
joint_axes(h,k);
