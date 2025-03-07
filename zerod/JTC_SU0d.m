% script de run de JTC_SU
scenar ='';
while isempty(scenar)
	scenar = upper(input('Scenario ? I(p High), H(ybrid), S(teady-state) or D(ay long) [ R(efrence 17266)]','s'));
end
scenar = upper(scenar(1));
switch scenar
case 'I'
        temps    = linspace(0,100,301)';
	b0    = 2.68;
	R     = 3.06;
	a     = 1.15;
	K95     = 1.76;
	d95     = 0.45;
	ip    = 5.5e6;
	nbar  = 0.4 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 0;
	picrh = 0;
	pecrh = 7e6;
	pnbi  = 34e6;
	zeff  = 2;
	li    = 0.7;
	hmore = 1.31;
	xece  = 0.6;
	rtang =R;
	einj  = 85e3;
	pnbi_p  = 24e6;
	pnbi_n  = 10e6;
	scenario = 'High Ip';
	ane = 3;
	vane = 1.01;
	xiioxie = 1;
    fpped = 1;	   
    angle_nbi = 55;	

case 'H'
        temps    = linspace(0,100,301)';
	b0    = 2.68;
	R     = 3.06;
	a     = 1.15;
	K95     = 1.76;
	d95     = 0.45;
	ip    = 3.5e6;
	nbar  = 0.85 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 0;
	picrh = 0;
	pecrh = 7e6;
	pnbi  = 34e6;
	zeff  = 2;
	li    = 0.7;
	hmore = 1.1;
	xece  = 0.6;
	rtang = R+a/4;
	einj  = 85e3;
	pnbi_p  = 24e6;
	pnbi_n  = 10e6;
	scenario = 'Hybrid';
	ane = 3;
	vane = 1.01;
	xiioxie = 1;
    angle_nbi = 55;	
    fpped = 1;	

case 'S'
        temps    = linspace(0,100,301)';
	b0    = 1.8;
	R     = 3.06;
	a     = 1.15;
	K95     = 1.76;
	d95     = 0.45;
	ip    = 2.4e6;
	nbar  = 0.86 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 0;
	picrh = 0;
	pecrh = 7e6;
	pnbi  = 34e6;
	zeff  = 2;
	li    = 0.7;
	hmore = 1.33;
	xece  = 0.6;
	rtang = R+a/4;
	einj  = 85e3;
	pnbi_p  = 24e6;
	pnbi_n  = 10e6;
	scenario = 'Steady state';
	ane = 3;
	vane = 1.01;
	xiioxie = 1;
    fpped = 1;	
    angle_nbi = 55;	

case 'D'
        temps    = logspace(0,4.46,301)';
	b0    = 1.3;
	R     = 3.06;
	a     = 1.15;
	K95     = 1.76;
	d95     = 0.45;
	ip    = 1e6;
	nbar  = 0.86 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 0;
	picrh = 0;
	pecrh = 7e6;
	pnbi  = 8e6;
	zeff  = 2;
	li    = 0.7;
	hmore = 1.31;
	xece  = 0.6;
	rtang = R+a/4;
	einj  = 85e3;
	pnbi_p  = 8e6;
	pnbi_n  = 0e6;
	scenario = 'Day long';
	ane = 3;
	vane = 1.01;
    xiioxie = 1;
    fpped = 1;	
     angle_nbi = 55;	

case 'R'
    
        temps    = linspace(0,100,301)';
	b0    = 4.4;
	R     = 3.05;
	a     = 0.71;
	K95     = 1.7;
	d95     = 0.2;
	ip    = 2.17e6;
	nbar  = 3.8e19;
	plh   = 0;
	picrh = 0;
	pecrh = 0e6;
	pnbi  = 33.3e6;
	zeff  = 2;
	li    = 0.7;
	hmore = 1.4;
	xece  = 0.6;
	rtang = 0;
	einj  = 90e3;
	pnbi_p  = 33e6;
	pnbi_n  = 0e6;
	scenario = 'reference: 17266';
	ane = 4;
	vane = 2.5;
    xiioxie = 0.9;    
    fpped = 0.3;
    angle_nbi = 20;	

end

% formule a partir de k0 = (K+1)/2 et K(x) = k0 + (K-k0) x^ 4
K  = 2 .* (K95 - 0.5 .* (1 - 0.95 .^4)) ./ (1 + 0.95 .^4);

d  = d95 ./ (0.95 .^ 2);

% limite de stabilite verticale
k95lim = 2.22 - 0.17 .* R ./ a;
if K95 > (1.05.*k95lim)
	fprintf('Warning : K too high, plamsa could be verticaly instable  (K95lim = %g  & K95 = %g)\n',k95lim,K95);
end

z0dinput = zerod_scalaire(temps,b0,R,a,K,d,ip,nbar,plh,picrh,pecrh,pnbi,zeff,xece,hmore,0,0,0, ...
                          'gaz',2,'frhe0',0,'tauhemul',5,'ane',ane,'vane',vane, ...
			  'scaling',0,'l2hscalin',0,'modeh',1,'l2hmul',0,'fpped',fpped, ...
			  'xiioxie',xiioxie,'kishape',3,'qdds',0.95,'kidds',3,'vloop',0,'runaway',1,'modeboot',1,'vref',0, ...
			  'zeff',4,'zmax',10,'zimp',6,'rimp',0.1,'matthews',1, ...
			  'frad',1,'rw',0.7,'angle_ece',90,'synergie',1, ...
			  'sens',1,'angle_nbi',angle_nbi,'einj',einj,'rtang',rtang,'lhmode',3,'etalh',0.8, ...
			  'npar0',2,'freqlh',5,'wlh',0,'xlh',0.5,'dlh',0.5,'fwcd',1, ...
			  'mino','H','cmin',0.1,'nphi',25,'freq',72,'sitb',0,'tae',0, ...
			  'smhd',0,'tmhd',inf,'rip',0,'fprad',1/3,'li',li);

z0dinput.cons.nbar = z0dinput.cons.nbar .* min(1, 0.1 + temps ./ 10);
z0dinput.cons.pnbi = z0dinput.cons.pnbi .* min(1, 0.1 + temps ./ 30);
z0dinput.cons.pecrh = z0dinput.cons.pecrh .* min(1, 0.3 + temps ./ 5);
z0dinput.cons.ip = z0dinput.cons.ip .* min(1, 0.1 + temps ./ 5);
z0dinput.cons.hmore = z0dinput.cons.hmore .* min(1, 0.1 + temps ./ 60);

%z0dinput.cons.iso = z0dinput.cons.iso .* min(1, 0.5 + temps ./ 500);


[zs,infovoid,profli] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);

post.z0dinput = z0dinput;
post.zerod    = zs;
post.profil0d =profli;
op0d =post.z0dinput.option;

rfan = (3.56e6 + 14.03e6) ./ 3.56e6 ;
padd = z0dinput.cons.picrh + z0dinput.cons.pecrh + z0dinput.cons.pnbi + zs.pohm + zs.plh;
betanth = zs.wth .* (1.6.*pi./3) .* z0dinput.geo.a ./ zs.vp ./ z0dinput.geo.b0 ./ zs.ip;
fprintf('b0 = %g, R = %g, a = %g, Ka = %g, da = %g, epsi = %g\n', z0dinput.geo.b0(end-2),  z0dinput.geo.R(end-2), ...
        z0dinput.geo.a(end-2),z0dinput.geo.K(end-2),z0dinput.geo.d(end-2),z0dinput.geo.R(end-2)./z0dinput.geo.a(end-2));
fprintf('Ip = %g, qa = %g, q95 = %g, q0 = %g, qmin = %g, li = %g\n', ...
        zs.ip(end-2)./1e6,zs.qa(end-2),zs.q95(end-2), ...
        zs.q0(end-2),zs.qmin(end-2),zs.li(end-2));
fprintf('nbar/ngr = %g,<nD> = %g, <nT> = %g, <nHe> =%g,<nimp> =%g\n',zs.nbar(end-2) ./1e20./  (zs.ip(end-2) / 1e6) .*  ...
        (pi.* z0dinput.geo.a(end-2) .^ 2),zs.nDm(end-2)./1e19,zs.nTm(end-2)./1e19,zs.nhem(end-2)./1e19,zs.nimpm(end-2))./1e19;
fprintf('taue = %g, tauhe = %g , HH = %g, Wth = %g, W = %g, Wfast = %g\n', ...
        zs.taue(end-2),zs.tauhe(end-2),zs.taue(end-2)./zs.tauh(end-2).*zs.hitb(end-2),  ...
	zs.wth(end-2)./1e6,zs.w(end-2)./1e6,(zs.w(end-2) -  zs.wth(end-2)) ./1e6);
fprintf('Plh = %g, Picrh = %g, Pecrh = %g, Pnbi = %g, Pline = %g, Pbrem = %g, Pcyclo = %g\n',zs.plh(end-2)/1e6, ...
        zs.picrh(end-2)/1e6,zs.pecrh(end-2)/1e6,zs.pnbi(end-2)/1e6,(zs.prad(end-2)+ zs.pradsol(end-2))/1e6,zs.pbrem(end-2)/1e6, zs.pcyclo(end-2) ./1e6);
fprintf('<ne> = %g, ne0 = %g, ne0/<ne> = %g, <Te> = %g , Te0 = %g, Te0/<Te> = %g,  Ti/Te = %g\n',zs.nem(end-2)/1e19,zs.ne0(end-2)/1e19, ...
         zs.ne0(end-2) ./ zs.nem(end-2) ,zs.tem(end-2)/1e3,zs.te0(end-2)/1e3,zs.te0(end-2) ./ zs.tem(end-2),zs.tite(end -2));
fprintf('Pfus_in = %g, Pfus = %g, Q = %g, betan(total) = %g , betan(th) = %g, frad = %g  \n',zs.pfus(end-2)/1e6,zs.pfus(end-2).* rfan./1e6, ...
         rfan .* zs.pfus(end-2) ./ padd(end-2),zs.betan(end-2).*100,betanth(end-2).*100, ...
	 (zs.prad(end-2) + zs.pbrem(end-2) + zs.pcyclo(end-2)) ./ zs.pin(end-2));
fprintf('Iboot/Ip = %g, Icd/Ip = %g, Vloop = %g , Zeff(line) = %g, li =%g\n',zs.iboot(end-2) ./zs.ip(end-2),zs.icd(end-2) ./ zs.ip(end-2), ...
        zs.vloop(end-2),zs.zeff(end-2)./zs.zmszl(end-2),zs.li(end-2));
fprintf('Vol = %g, Section = %g, Sext = %g , Llcms = %g\n',zs.vp(end-2),zs.sp(end-2),zs.sext(end-2),zs.peri(end-2));
fprintf('frHe = %g, frImp = %g @ Z= %d  &  %g @ Z = %d, <Zeff> = %g \n', ...
         zs.nhem(end-2)./zs.nem(end-2),zs.nimpm(end-2)./zs.nem(end-2),op0d.zimp, ...
	 zs.nimpm(end-2)./zs.nem(end-2).* op0d.rimp,op0d.zmax, ...
	 (zs.n1m(end-2) + 4 .* zs.nhem(end-2) +zs.nimpm(end-2) .* op0d.zimp .^ 2 + ...
	 zs.nimpm(end-2) .* op0d.zmax  .^ 2 .* op0d.rimp) ./ zs.nem(end-2) );


z0plotsc

% calcul des neutrons par type d'injecteur
option.einj = 85e3;
[void1,void2,zs.ndd_nbi_thp] = z0neutron_dd(option,cons,zs,profli,pnbi_p); 
option.einj = 500e3;
[void1,void2,zs.ndd_nbi_thn] = z0neutron_dd(option,cons,zs,profli,pnbi_n);

% correction du beam-beam
fprintf('BB pp:');
zs.ndd_nbi_nbi_p  = z0beambeam(zs.temps,profli.xli,profli.vpr,profli.nep,profli.tep,profli.n1p,profli.tip,profli.zeff,profli.pnbi, ....
                            zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,zs.nem,85e3,85e3,pnbi_p,pnbi_p,zs.mu0_nbi,1);
fprintf('BB nn:');
zs.ndd_nbi_nbi_n  = z0beambeam(zs.temps,profli.xli,profli.vpr,profli.nep,profli.tep,profli.n1p,profli.tip,profli.zeff,profli.pnbi, ....
                            zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,zs.nem,500e3,500e3,pnbi_n,pnbi_n,zs.mu0_nbi,1);
fprintf('BB pn:');
zs.ndd_nbi_nbi_pn  = z0beambeam(zs.temps,profli.xli,profli.vpr,profli.nep,profli.tep,profli.n1p,profli.tip,profli.zeff,profli.pnbi, ....
                            zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,zs.nem,85e3,500e3,pnbi_p,pnbi_n,zs.mu0_nbi,0);
%  fprintf('BB total:');
%  zs.ndd_nbi_nbi_tot  = z0beambeam_test(zs.temps,profli.xli,profli.vpr,profli.nep,profli.tep,profli.n1p,profli.tip,profli.zeff,profli.pnbi, ....
%                              zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,zs.nem,85e3,500e3,pnbi_p,pnbi_n,zs.mu0_nbi,1);



h = findobj(0,'type','figure','tag','neutron');
if isempty(h)
       h=figure('tag','neutron');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

zs.ndd_total = zs.ndd_th + zs.ndd_nbi_thp + zs.ndd_nbi_thn + zs.ndd_nbi_nbi_p + zs.ndd_nbi_nbi_n + zs.ndd_nbi_nbi_pn;
               
 
plot(zs.temps,zs.ndd_total,'k', ... 
     zs.temps,zs.ndd_th,'g', ...
      zs.temps,zs.ndd_nbi_thp,'b',zs.temps,zs.ndd_nbi_thn,'c', ...
      zs.temps,zs.ndd_nbi_nbi_p,'r',zs.temps,zs.ndd_nbi_nbi_n,'m',zs.temps,zs.ndd_nbi_nbi_pn,'k:')
      
xlabel('time (s)');
ylabel('neutron/s');
title(sprintf('JT60 SA (%s)',scenario));
legend('Total','Thermal','Posifif : beam-plasma','Negative: beam-plasma','Posifif : beam-beam','Negative: beam-beam','Negative-Positive: beam-beam')
      
    
