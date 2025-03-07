% plasma parameters
%  R   = 1.8;
%  a   = 0.6;
%  K95     = 1.7;
%  d95     = 0.4;
%  b0  = 6.7;
%  ip  = 6e6;
%  picrh = 30e6;
%  pnbi  = 0;
%  pecrh = 4e6;
%  plh   = 6e6;
%  zeff  = 1.3;
%  li    = 0.7;
%  hmore = 1;
%  xece  = 0.6;
%  xiioxie = 0;
%  fpped = 1;	   




scenar ='';
while isempty(scenar) 
	scenar = upper(input('Scenario ? H1, H2, HYB1 (HYB11) , HYB2 (HYB12), SS1, SS2 or SS3 ?','s'));
end
switch scenar
case 'H1'
        tmax    = 12;
	ipmax   = 5e6;
	nbar    = 1.7e20;
	R       = 1.8;
	a       = 0.6;
	K95     = 1.7;
	d95     = 0.4;
	b0      = 6.7;
	picrh   = 20e6;
	pnbi    = 0;
	pecrh   = 0e6;
	plh     = 0e6;
	zeff    = 1.2;
	li      = 0.7;
	hmore   = 1;
	xece    = 0.6;
	xiioxie = 0;
	fpped   = 0.5;	
	sitb    = 0;
	   
case 'H2'
        tmax    = 12;
	ipmax   = 5e6;
	nbar    = 2.6e20;
	R       = 1.8;
	a       = 0.6;
	K95     = 1.7;
	d95     = 0.4;
	b0      = 6.7;
	picrh   = 20e6;
	pnbi    = 0;
	pecrh   = 0e6;
	plh     = 0e6;
	zeff    = 1.18;
	li      = 0.7;
	hmore   = 1;
	xece    = 0.6;
	xiioxie = 0;
	sitb    = 0;
	fpped   = 0.5;	   
	
case 'HYB1'
        tmax    = 20;
	ipmax   = 3.6e6;
	nbar    = 1.95e20;
	R       = 1.98;
	a       = 0.76;
	K95     = 1.74;
	d95     = 0.53;
	b0      = 6.25;
	picrh   = 20e6;
	pnbi    = 0;
	pecrh   = 0e6;
	plh     = 0e6;
	zeff    = 1.2;
	li      = 0.7;
	hmore   = 1.3;
	xece    = 0.6;
	xiioxie = 0;
	fpped   = 0.5;	   
	sitb    = 0;
	
case 'HYB11'
        tmax    = 20;
	ipmax   = 3.6e6;
	nbar    = 1.3e20;
	R       = 1.98;
	a       = 0.76;
	K95     = 1.74;
	d95     = 0.53;
	b0      = 6.25;
	picrh   = 20e6;
	pnbi    = 0;
	pecrh   = 0e6;
	plh     = 0e6;
	zeff    = 1.2;
	li      = 0.7;
	hmore   = 1.3;
	xece    = 0.6;
	xiioxie = 0;
	sitb    = 0;
	fpped   = 0.5;	   
	
case 'HYB2'
        tmax    = 20;
	ipmax   = 3.6e6;
	nbar    = 1.3e20;
	R       = 1.98;
	a       = 0.76;
	K95     = 1.74;
	d95     = 0.53;
	b0      = 6.25;
	picrh   = 30e6;
	pnbi    = 0;
	pecrh   = 4e6;
	plh     = 6e6;
	zeff    = 1.2;
	li      = 0.7;
	hmore   = 1.3;
	xece    = 0.6;
	xiioxie = 0;
	fpped   = 0.5;	   
	sitb    = 0;
	
case 'HYB12'
        tmax    = 20;
	ipmax   = 3.6e6;
	nbar    = 1.95e20;
	R       = 1.98;
	a       = 0.76;
	K95     = 1.74;
	d95     = 0.53;
	b0      = 6.25;
	picrh   = 30e6;
	pnbi    = 0;
	pecrh   = 4e6;
	plh     = 6e6;
	zeff    = 1.2;
	li      = 0.7;
	hmore   = 1.3;
	xece    = 0.6;
	xiioxie = 0;
	fpped   = 0.5;	   
	sitb    = 0;
	
case 'SS1'
        tmax    = 200;
	ipmax   = 2.8e6;
	nbar    = 0.9e20;
	R       = 2;
	a       = 0.74;
	K95     = 1.56;
	d95     = 0.49;
	b0      = 6;
	picrh   = 20e6;
	pnbi    = 0;
	pecrh   = 4e6;
	plh     = 6e6;
	zeff    = 1.46;
	li      = 0.7;
	hmore   = 1.5;
	xece    = 0.3;
	xiioxie = 0;
	fpped   = 0.5;	   
	sitb    = 0;
	
case 'SS2'
        tmax    = 200;
	ipmax   = 2.8e6;
	nbar    = 1.35e20;
	R       = 2;
	a       = 0.74;
	K95     = 1.56;
	d95     = 0.49;
	b0      = 6;
	picrh   = 20e6;
	pnbi    = 0;
	pecrh   = 4e6;
	plh     = 6e6;
	zeff    = 1.46;
	li      = 0.7;
	hmore   = 1.5;
	xece    = 0.3;
	xiioxie = 0;
	fpped   = 0.5;	   
	sitb    = 0;
	
case 'SS3'
        tmax    = 200;
	ipmax   = 2.8e6;
	nbar    = 1.9e20;
	R       = 2;
	a       = 0.74;
	K95     = 1.56;
	d95     = 0.49;
	b0      = 6;
	picrh   = 20e6;
	pnbi    = 0;
	pecrh   = 4e6;
	plh     = 6e6;
	zeff    = 1.46;
	li      = 0.7;
	hmore   = 1.5;
	xece    = 0.3;
	xiioxie = 0;
	fpped   = 0.5;	  
	sitb    = 0;
	 
case 'SS11'
        tmax    = 200;
	ipmax   = 2.8e6;
	nbar    = 0.9e20;
	R       = 2;
	a       = 0.74;
	K95     = 1.56;
	d95     = 0.49;
	b0      = 6;
	picrh   = 20e6;
	pnbi    = 0;
	pecrh   = 4e6;
	plh     = 6e6;
	zeff    = 1.46;
	li      = 0.7;
	hmore   = 1;
	xece    = 0.3;
	xiioxie = 0;
	fpped   = 0.5;	   
	sitb    = 3;
	
case 'SS12'
        tmax    = 200;
	ipmax   = 2.8e6;
	nbar    = 1.35e20;
	R       = 2;
	a       = 0.74;
	K95     = 1.56;
	d95     = 0.49;
	b0      = 6;
	picrh   = 20e6;
	pnbi    = 0;
	pecrh   = 4e6;
	plh     = 6e6;
	zeff    = 1.46;
	li      = 0.7;
	hmore   = 1;
	xece    = 0.3;
	xiioxie = 0;
	fpped   = 0.5;	   
	sitb    = 3;
	
case 'SS13'
        tmax    = 200;
	ipmax   = 2.8e6;
	nbar    = 1.9e20;
	R       = 2;
	a       = 0.74;
	K95     = 1.56;
	d95     = 0.49;
	b0      = 6;
	picrh   = 20e6;
	pnbi    = 0;
	pecrh   = 4e6;
	plh     = 6e6;
	zeff    = 1.46;
	li      = 0.7;
	hmore   = 1;
	xece    = 0.3;
	xiioxie = 0;
	fpped   = 0.5;	  
	sitb    = 3;
	 
otherwise
	error('Unknown scenario')
end

% script pour le test de FT3 
% formule a partir de k0 = (K+1)/2 et K(x) = k0 + (K-k0) x^ 4
K  = 2 .* (K95 - 0.5 .* (1 - 0.95 .^4)) ./ (1 + 0.95 .^4);
d  = d95 ./ (0.95 .^ 2);

% limite de stabilite verticale
k95lim = 2.22 - 0.17 .* R ./ a;
if K95 > (1.05.*k95lim)
	fprintf('Warning : K too high, plamsa could be verticaly instable  (K95lim = %g  & K95 = %g)\n',k95lim,K95);
end

temps = linspace(0.5,tmax,301)';
z0dinput = zerod_scalaire(temps,b0,R,a,K,d,ipmax,nbar,plh,picrh,pecrh,pnbi,zeff,xece,hmore,0,0,0, ...
                          'gaz',2,'frhe0',0,'tauhemul',5,'ane',0,'vane',1.1, ...
			  'scaling',0,'l2hscalin',0,'modeh',1,'l2hmul',0,'fpped',fpped, ...
			  'xiioxie',xiioxie,'kishape',3,'qdds',0.8,'kidds',3,'vloop',0,'runaway',1,'modeboot',1,'vref',0, ...
			  'zeff',0,'zmax',74,'zimp',18,'rimp',0.1,'matthews',0, ...
			  'frad',1,'rw',0.7,'angle_ece',90,'synergie',1, ...
			  'sens',1,'angle_nbi',60,'einj',1e5,'rtang',R+a/4,'lhmode',0,'etalh',0.8, ...
			  'npar0',2,'freqlh',5,'wlh',0.3,'xlh',0.5,'dlh',0.5,'fwcd',0, ...
			  'mino','He3','cmin',0.02,'nphi',8,'freq',68,'sitb',sitb,'tae',0, ...
			  'smhd',0,'tmhd',0,'rip',0,'fprad',1/3,'li',li);

			  
% evolution temporelle	
tv     = [0.5	1.5	2.5	4.5		6	tmax];
ipv    = [0.5e6	2e6	3e6	ipmax.*4./5	ipmax	ipmax];
nv     = [4e19	4e19	max(4e19,nbar/3)	nbar.*4./5	nbar    nbar];
picrhv = [0     0       0       0.5             1       1   ] .* picrh;
plhv   = [0     0       0.5       1               1       1   ] .* plh;
pecrhv = [0     1       1       1               1       1   ] .* pecrh;
pnbiv  = [0     0       0       0               1       0   ] .* pnbi;
Rv     = [1     1       1       1               1       1   ] .* R;
av     = [0.5   1       1       1               1       1   ] .* a;
Kv     = [0     0       0.5     1               1       1   ] .* (K - 1) + 1;
dv     = [0     0       0       1               1       1   ] .* d;
		  
z0dinput.cons.ip = pchip(tv,ipv,temps);
z0dinput.cons.picrh = pchip(tv,picrhv,temps);
z0dinput.cons.pnbi = pchip(tv,pnbiv,temps);
z0dinput.cons.plh = pchip(tv,plhv,temps);
z0dinput.cons.pecrh = pchip(tv,pecrhv,temps);
z0dinput.geo.a   = pchip(tv,av,temps);
z0dinput.geo.R   = pchip(tv,Rv,temps);
z0dinput.geo.K   = pchip(tv,Kv,temps);
z0dinput.geo.d   = pchip(tv,dv,temps);
z0dinput.cons.nbar = pchip(tv,nv,temps);

%z0dinput.cons.iso = z0dinput.cons.iso .* min(1, 0.5 + temps ./ 500);
z0dinput.cons.zeff  = zeff  + 2.3 .*  (min(z0dinput.cons.nbar) ./ z0dinput.cons.nbar) .^ 2.6;
z0dinput.cons.hmore = (z0dinput.cons.hmore-1) .* min(1,max(0,1 + (temps - 4.5) ./ (7 - 4.5))) +1;

% utilisation des moments
if isfield(z0dinput.exp0d,'Rsepa')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Rsepa');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Zsepa');
end


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
save(sprintf('METIS_FT3_%s',scenar),'post');