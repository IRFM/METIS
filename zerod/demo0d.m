% script de run de demo
clear sepa_option

model ='';
while isempty(model)
	model = upper(input('model ? {A,B,C,D, I(ter), T(est), G(reat Iter), H(Iter Ohmique) , S(no CS)} ','s'));
end
mode ='';
while isempty(mode)
	mode = upper(input('mode de calcul ? {P(PCS) N(atural)} ','s'));
end
switch model
case 'I'
	b0    = 5.3;
	R     = 6.2;
	a     = 2;
	K95     = 1.7;
	d95     = 0.33;
	ip    = 15e6;
	nbar  = 0.85 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 0;
	picrh = 17e6;
	pecrh = 0e6;
	pnbi  = 33e6;
	zeff  = 1.5;
	li    = 0.7;
	hmore = 1.0;
	iso   = 1;
	ftnbi = 0;
	ane   = 1.01;
	einj  = 1e6;
	frad  = 0.6;
	rw    = 0.7;
	fpped = 1;
	fprad  = 1/3;

	sepa_option.rxup      = 0.466;     % upper triangularity (minor radius unit)
	sepa_option.zxup      = 1.687;    % upper altitude X point (minor radius unit)
	sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
	sepa_option.ra        = 6.2;       % major radius R0 (m) [6.2]
	sepa_option.za        = 0.65;       % altitude of the magnetic axis (m) [0.9]
	sepa_option.a         = 2;         % minor radius (m) [2]
	sepa_option.rxdo      = 0.568;     % lower triangularity (minor radius unit)
	sepa_option.zxdo      = 2.001;       % lower altitude X point (minor radius unit)
	sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
	sepa_option.b0        = 5.3 ./ (1 - 2 ./ 6.2 - 1 / 6.2);      % magnetic field at R0
	sepa_option.delta     = 1;      % magnetic field at R0
	sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
	sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]



case 'H'
	b0    = 5.3;
	R     = 6.2;
	a     = 2;
	K95     = 1.7;
	d95     = 0.33;
	ip    = 15e6;
	nbar  = 0.5 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 0;
	picrh = 0.*20e6;
	pecrh = 0e6;
	pnbi  = 0e6;
	zeff  = 1.23;
	li    = 0.7;
	hmore = 1.0;
	iso   = 1;
	ftnbi = 0.5;
	ane   = 1.01;
	einj  = 1e6;
	frad  = 1;
	rw    = 0.7;
	fpped = 1;
	fprad  = 1/3;

	sepa_option.rxup      = 0.466;     % upper triangularity (minor radius unit)
	sepa_option.zxup      = 1.687;    % upper altitude X point (minor radius unit)
	sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
	sepa_option.ra        = 6.2;       % major radius R0 (m) [6.2]
	sepa_option.za        = 0.65;       % altitude of the magnetic axis (m) [0.9]
	sepa_option.a         = 2;         % minor radius (m) [2]
	sepa_option.rxdo      = 0.568;     % lower triangularity (minor radius unit)
	sepa_option.zxdo      = 2.001;       % lower altitude X point (minor radius unit)
	sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
	sepa_option.b0        = 5.3 ./ (1 - 2 ./ 6.2 - 1 / 6.2);      % magnetic field at R0
	sepa_option.delta     = 1;      % magnetic field at R0
	sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
	sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]



case 'G'
	b0    = 5.68;
	R     = 8.14;
	a     = 2.80;
	K95   = 1.6;
	d95   = 0.24;
	ip    = 21e6;
	nbar  = 1.1 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 0;
	picrh = 0.*20e6;
	pecrh = 50e6;
	pnbi  = 100e6;
	zeff  = 1.24;
	li    = 0.7;
	hmore = 1.0;
	iso   = 1;
	ftnbi = 0.5;
	ane   = 1;
	einj  = 1e6;
	frad  = 0.6;
	rw    = 0.7;
	fpped = 1;
	fprad  = 1/3;
	sepa_option = [];

case 'A'
	b0    = 7;
	R     = 9.55;
	a     = 9.55 /3;
	K95   = 1.7 ;
	d95   = 0.25;
	ip    = 30.5e6;
	nbar  = 1.2 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 0;
	picrh = 176e6;
	pecrh = 0;
	pnbi  = 70e6;
	zeff  = 1.24;
	li    = 0.7;
	hmore = 1.2;
	iso   = 1;
	ftnbi = 0.5;
	ane   = 1.3;
	einj  = 2e6;
	frad   = 1;
	rw     = 0.7;
	fpped = 1;
	fprad  = 0;

	sepa_option.rxup      = 0.353;     % upper triangularity (minor radius unit)
	sepa_option.zxup      = 1.687;    % upper altitude X point (minor radius unit)
	sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
	sepa_option.ra        = 9.55;       % major radius R0 (m) [6.2]
	sepa_option.za        = 0.0;       % altitude of the magnetic axis (m) [0.9]
	sepa_option.a         = 9.55/3;         % minor radius (m) [2]
	sepa_option.rxdo      = 0.430;     % lower triangularity (minor radius unit)
	sepa_option.zxdo      = 2.001;       % lower altitude X point (minor radius unit)
	sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
	sepa_option.b0        = 13.1;      % magnetic field at R0
	sepa_option.delta     = 1.26;      % magnetic field at R0
	sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
	sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]



case 'B'
	b0    = 6.9;
	R     = 8.6;
	a     = 8.6 /3;
	K95     = 1.7;
	d95     = 0.25;
	ip    = 28e6;
	nbar  = 1.2 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 0;
	picrh = 200e6;
	pecrh = 0;
	pnbi  = 70e6;
	zeff  = 1.24;
	li    = 0.7;
	hmore = 1.2;
	iso   = 1;
	ftnbi = 0.5;
	ane   = 1.3;
	einj  = 1e6;
	frad   = 0.6;
	rw     = 0.7;
	fprad  = 0
	fpped = 1;

	sepa_option = [];

case 'C'
	b0    = 6;
	R     = 7.4;
	a     = R /3;
	K95     = 1.87;
	d95     = 0.47;
	ip    = 18.5e6;
	nbar  = 1.3 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 100e6;
	picrh = 0e6;
	pecrh = 0e6;
	pnbi  = 80e6;
	zeff  = 1.8;
	li    = 0.7;
	hmore = 1.37
	iso   = 1;
	ftnbi = 0.5;
	ane   = 1.1;
	einj  = 1e6;
	frad  = 1;
	rw    = 0.7;
	fprad  = 1/3;
	fpped = 1;

	sepa_option.rxup      = 0.6309;     % upper triangularity (minor radius unit)
	sepa_option.zxup      = 1.921;    % upper altitude X point (minor radius unit)
	sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
	sepa_option.ra        = 7.4;       % major radius R0 (m) [6.2]
	sepa_option.za        = 0;       % altitude of the magnetic axis (m) [0.9]
	sepa_option.a         = 7.4/3        % minor radius (m) [2]
	sepa_option.rxdo      = 0.7691;     % lower triangularity (minor radius unit)
	sepa_option.zxdo      = 2.279;       % lower altitude X point (minor radius unit)
	sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
	sepa_option.b0        = 13.6;      % magnetic field at R0
	sepa_option.delta     = 1.69;
	sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
	sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]

case 'D'
	b0    = 5.6;
	R     = 6.1;
	a     = 6.1 /3;
	K95     = 1.9;
	d95     = 0.47;
	ip    = 14.1e6;
	nbar  = 1.5 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 0e6;
	picrh = 0e6;
	pecrh = 38e6;
	pnbi  = 33e6;
	zeff  = 1.24;
	li    = 0.7;
	hmore = 1.5;
	iso   = 1;
	ftnbi = 0.5;
	ane   = 1.5;
	einj  = 1e6;
	frad  = 1;
	rw    = 0.7;
	fpped = 1;
	fprad  = 0;
	sepa_option = [];

case 'T'
	R = [];
	while(isempty(R))
		R     = input('R (m) ? ');
	end
	rsa = [];
	while(isempty(rsa))
		rsa     = input('R/a (su) ? ');
	end
	a     = R / rsa;
	b0    = (5.3 .* 6.2) ./(6.2-2-1) .*(R-a-1) ./ R;
	K95   = 2.22 - 0.17 .* R ./ a;
	d95   = 0.33;
	ip    = 5e6 .* a .^ 2 .* b0 .* (1+ K95 .^ 2)./2 ./ R ./3.* (1+2 .* (a./R).^2) .*...
	        (1.24 - 0.54 .* K95 + 0.3 .* (K95 .^ 2 + d95 .^ 2) + 0.13 .* d95);
	nbar  = 0.95 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 0e6;
	picrh = 0e6;
	pecrh = 0;
	pnbi  = 0.042e6 .* (nbar/1e20) .^ 0.73 .* b0 .^0.74 .* (K95.*a.*R.*4.*pi.^2) .^0.98./2;
	zeff  = 1.24;
	li    = 0.7;
	hmore = 1;
	iso   = 1;
	ftnbi = 0.5;
	ane   = 1;
	einj  = 1e6;
	frad  = 0.6;
	rw    = 0.7;
	fpped = 0.8;
	sepa_option = [];
	fprad  = 1/3;
case 'S'
	b0    = 6;
	R     = 5.5;
	a     = 2.1;
	K95     = 2;
	d95     = 0.4;
	ip    = 16.7e6;
	nbar  = 1.15 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
	plh   = 000e6;
	picrh = 0e6;
	pecrh = 80e6;
	pnbi  = 20e6;
	zeff  = 1.5;
	li    = 0.1;
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


end

% formule a partir de k0 = (K+1)/2 et K(x) = k0 + (K-k0) x^ 4
K  = 2 .* (K95 - 0.5 .* (1 - 0.95 .^4)) ./ (1 + 0.95 .^4);

d  = d95 ./ (0.95 .^ 2);

% limite de stabilite verticale
k95lim = 2.22 - 0.17 .* R ./ a;
if K95 > (1.05.*k95lim)
	fprintf('Warning : K too high, plamsa could be verticaly instable  (K95lim = %g  & K95 = %g)\n',k95lim,K95);
end

switch model
case 'H'
	temps    = linspace(1,301,301)';
otherwise
	temps    = logspace(0,4,301)';
end

xece = 0.3;

switch mode
case   'P'
	z0dinput = zerod_scalaire(temps,b0,R,a,K,d,ip,nbar,plh,picrh,pecrh,pnbi,zeff,xece,hmore,iso,ftnbi,0, ...
                          'gaz',3,'frhe0',0,'tauhemul',5,'ane',4,'vane',ane, ...
			  'scaling',0,'l2hscalin',0,'modeh',1,'l2hmul',0,'fpped',fpped, ...
			  'xiioxie',2,'kishape',3,'qdds',0.95,'kidds',3,'vloop',0,'runaway',0,'modeboot',1,'vref',0, ...
			  'zeff',0,'zmax',18,'zimp',4,'rimp',0.06,'matthews',0, ...
			  'frad',frad,'rw',rw,'angle_ece',90,'synergie',1, ...
			  'sens',1,'angle_nbi',50,'einj',einj,'rtang',R+a/4,'lhmode',3,'etalh',0.8, ...
			  'npar0',2,'freqlh',5,'wlh',a/2,'xlh',0,'dlh',0,'fwcd',0, ...
			  'mino','T','cmin',1,'nphi',25,'freq',72,'sitb',0,'tae',0, ...
			  'smhd',0,'tmhd',inf,'rip',0,'fprad',fprad,'li',1,'configuration',2);
case 'N'
	z0dinput = zerod_scalaire(temps,b0,R,a,K,d,ip,nbar,plh,picrh,pecrh,pnbi,zeff,xece,ones(size(hmore)),iso,ftnbi,0, ...
                          'gaz',3,'frhe0',0,'tauhemul',5,'ane',3,'vane',1, ...
			  'scaling',0,'l2hscalin',0,'modeh',1,'l2hmul',0,'fpped',1, ...
			  'xiioxie',0,'kishape',3,'qdds',0.95,'kidds',3,'vloop',0,'runaway',1, ...
			  'modeboot',1,'vref',0, ...
			  'zeff',2,'zmax',10,'zimp',4,'rimp',0.06,'matthews',0, ...
			  'frad',1,'rw',rw,'angle_ece',90,'synergie',0, ...
			  'sens',1,'angle_nbi',50,'einj',einj,'rtang',R,'lhmode',3,'etalh',0.8, ...
			  'npar0',2,'freqlh',5,'wlh',a/2,'xlh',0,'dlh',0,'fwcd',1, ...
			  'mino','T','cmin',1,'nphi',25,'freq',72,'sitb',3,'tae',1, ...
			  'smhd',0,'tmhd',9500,'rip',0,'fprad',fprad,'li',1,'configuration',2);

end
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
