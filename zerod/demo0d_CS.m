% script de run de demo
	b0    = 6;
	R     = 5.5;
	a     = 2.1;
	K95     = 2;
	d95     = 0.4;
	ip    = 16.7e6;
	nbar  = 1.15 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 000e6;
	picrh = 0e6;
	pecrh = 50e6;
	pnbi  = 50e6;
	zeff  = 1.5;
	li    = 0.5;
	hmore = 1.3;
	iso   = 1;
	ftnbi = 0;
	ane   = 1.6;
	einj  = 0.5e6;
	frad  = 1;
	rw    = 0.7;
	fprad  = 1/3;
	fpped = 1;
	sepa_option = [];
	xece = 0.3;



% formule a partir de k0 = (K+1)/2 et K(x) = k0 + (K-k0) x^ 4
K  = 2 .* (K95 - 0.5 .* (1 - 0.95 .^4)) ./ (1 + 0.95 .^4);

d  = d95 ./ (0.95 .^ 2);

% limite de stabilite verticale
k95lim = 2.22 - 0.17 .* R ./ a;
if K95 > (1.05.*k95lim)
	fprintf('Warning : K too high, plamsa could be verticaly instable  (K95lim = %g  & K95 = %g)\n',k95lim,K95);
end

temps    = logspace(0,4,301)';

xece = 0.3;

	z0dinput = zerod_scalaire(temps,b0,R,a,K,d,ip,nbar,plh,picrh,pecrh,pnbi,zeff,xece,hmore,iso,ftnbi,0, ...
                          'gaz',3,'frhe0',0,'tauhemul',5,'ane',5,'vane',1.2, ...
			  'scaling',0,'l2hscalin',0,'modeh',1,'l2hmul',0,'fpped',fpped, ...
			  'xiioxie',0,'kishape',3,'qdds',0.95,'kidds',3,'vloop',0,'runaway',0,'modeboot',1,'vref',0, ...
			  'zeff',1,'zmax',18,'zimp',4,'rimp',0.06,'matthews',0, ...
			  'frad',frad,'rw',rw,'angle_ece',90,'synergie',1, ...
			  'sens',1,'angle_nbi',50,'einj',einj,'rtang',R+a/4,'lhmode',3,'etalh',0.8, ...
			  'npar0',2,'freqlh',5,'wlh',a/2,'xlh',0,'dlh',0,'fwcd',0, ...
			  'mino','T','cmin',1,'nphi',25,'freq',72,'sitb',3,'tae',1, ...
			  'smhd',0,'tmhd',0,'rip',0,'fprad',fprad,'li',1,'configuration',2);

z0dinput.cons.nbar = z0dinput.cons.nbar .* min(1, 0.1 + temps ./ 300);
%z0dinput.cons.iso = z0dinput.cons.iso .* min(1, 0.5 + temps ./ 500);

if ~isempty(sepa_option)
	z0dinput = z0separatrix(z0dinput,sepa_option);
elseif isfield(z0dinput.exp0d,'Rsepa')
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Rsepa');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Zsepa');
end


[zs,infovoid,profli] = zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);

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
        zs.picrh(end-2)/1e6,zs.pecrh(end-2)/1e6,zs.pnbi(end-2)/1e6,(zs.prad(end-2) +zs.pradsol(end-2)) /1e6,zs.pbrem(end-2)/1e6, zs.pcyclo(end-2) ./1e6);
fprintf('<ne> = %g, ne0 = %g, ne0/<ne> = %g, <Te> = %g , Te0 = %g, Te0/<Te> = %g,  Ti/Te = %g\n',zs.nem(end-2)/1e19,zs.ne0(end-2)/1e19, ...
         zs.ne0(end-2) ./ zs.nem(end-2) ,zs.tem(end-2)/1e3,zs.te0(end-2)/1e3,zs.te0(end-2) ./ zs.tem(end-2),zs.tite(end -2));
fprintf('Palpha = %g, Pfus = %g, Q = %g, betan(total) = %g , betan(th) = %g, frad = %g  \n',zs.pfus(end-2)/1e6,zs.pfus(end-2).* rfan./1e6, ...
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
