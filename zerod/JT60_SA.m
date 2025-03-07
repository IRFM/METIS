% script de run de JTC_SA
% choix du gaz
gaz_liste ={'H','D','DT','He4'};
gaz = menu('JT60_SA gas',gaz_liste);

% choix de la forme du plasma en plateau
shape_liste={'Iter like CDR','Demo like CDR','Iter like reduce field','Demo like reduce field'};
shape = menu('JT60_SA plasma shape',shape_liste);
switch shape
case 1
	sepa_option.rxup      = 0.466;     % upper triangularity (minor radius unit)
	sepa_option.zxup      = 1.687;    % upper altitude X point (minor radius unit)
	sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
	sepa_option.ra        = 3.05;       % major radius R0 (m) [6.2]
	sepa_option.za        = 0;       % altitude of the magnetic axis (m) [0.9]
	sepa_option.a         = 0.984;         % minor radius (m) [2]
	sepa_option.rxdo      = 0.568;     % lower triangularity (minor radius unit)
	sepa_option.zxdo      = 2.001;       % lower altitude X point (minor radius unit)
	sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
	sepa_option.delta     = 1;      % magnetic field at R0
	sepa_option.b0        = 2.69 ./ (1 - sepa_option.a ./ sepa_option.ra - sepa_option.delta / sepa_option.ra);      % magnetic field at R0
	sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
	sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
	
	mode = 'cdr';

case 2
	sepa_option.rxup      = 0.466;     % upper triangularity (minor radius unit)
	sepa_option.zxup      = 1.72;    % upper altitude X point (minor radius unit)
	sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
	sepa_option.ra        = 3.06;       % major radius R0 (m) [6.2]
	sepa_option.za        = 0;       % altitude of the magnetic axis (m) [0.9]
	sepa_option.a         = 1.15;         % minor radius (m) [2]
	sepa_option.rxdo      = 0.568;     % lower triangularity (minor radius unit)
	sepa_option.zxdo      = 2.04;       % lower altitude X point (minor radius unit)
	sepa_option.apdo      = 19;   % lower separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amdo      = 71;   % lower separatrix angle (-R,X)  (HFS, degrees)
	sepa_option.delta     = 1;      % magnetic field at R0
	sepa_option.b0        = 2.27 ./ (1 - sepa_option.a ./ sepa_option.ra - sepa_option.delta / sepa_option.ra);      % magnetic field at R0
	sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
	sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]

	mode = 'cdr';

case 3
	sepa_option.rxup      = 0.466;     % upper triangularity (minor radius unit)
	sepa_option.zxup      = 1.687;    % upper altitude X point (minor radius unit)
	sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
	sepa_option.ra        = 2.94;       % major radius R0 (m) [6.2]
	sepa_option.za        = 0;       % altitude of the magnetic axis (m) [0.9]
	sepa_option.a         = 0.948;         % minor radius (m) [2]
	sepa_option.rxdo      = 0.568;     % lower triangularity (minor radius unit)
	sepa_option.zxdo      = 2.001;       % lower altitude X point (minor radius unit)
	sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
	sepa_option.delta     = 1;      % magnetic field at R0
	sepa_option.b0        = 2.68 ./ (1 - sepa_option.a ./ sepa_option.ra - sepa_option.delta / sepa_option.ra);      % magnetic field at R0
	sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
	sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
	mode = 'reduce';
	

case 4
	sepa_option.rxup      = 0.466;     % upper triangularity (minor radius unit)
	sepa_option.zxup      = 1.77;    % upper altitude X point (minor radius unit)
	sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
	sepa_option.ra        = 2.95;       % major radius R0 (m) [6.2]
	sepa_option.za        = 0;       % altitude of the magnetic axis (m) [0.9]
	sepa_option.a         = 1.18;         % minor radius (m) [2]
	sepa_option.rxdo      = 0.568;     % lower triangularity (minor radius unit)
	sepa_option.zxdo      = 2.11;       % lower altitude X point (minor radius unit)
	sepa_option.apdo      = 23;   % lower separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amdo      = 67;   % lower separatrix angle (-R,X)  (HFS, degrees)
	sepa_option.delta     = 1;      % magnetic field at R0
	sepa_option.b0        = 2.26 ./ (1 - sepa_option.a ./ sepa_option.ra - sepa_option.delta / sepa_option.ra);      % magnetic field at R0
	sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
	sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
	
	mode = 'reduce';

case 5
	sepa_option.rxup      = 0.5928;     % upper triangularity (minor radius unit)
	sepa_option.zxup      = 2.014;    % upper altitude X point (minor radius unit)
	sepa_option.apup      = 27.792;       % upper separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amup      = 64.13;       % upper separatrix angle (-R,X) (HFS, degrees)
	sepa_option.ra        = 2.941;       % major radius R0 (m) [6.2]
	sepa_option.za        = 0;       % altitude of the magnetic axis (m) [0.9]
	sepa_option.a         = 1.169;         % minor radius (m) [2]
	sepa_option.rxdo      = 0.5046;     % lower triangularity (minor radius unit)
	sepa_option.zxdo      = 1.963;       % lower altitude X point (minor radius unit)
	sepa_option.apdo      = 28.273;   % lower separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amdo      = 61.180;   % lower separatrix angle (-R,X)  (HFS, degrees)
	sepa_option.delta     = 1;      % magnetic field at R0
	sepa_option.b0        = 2.26 ./ (1 - sepa_option.a ./ sepa_option.ra - sepa_option.delta / sepa_option.ra);      % magnetic field at R0
	sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
	sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
	
	mode = 'reduce';

end

% scenario cible
scenario_liste ={'H mode @ fgr = 0.6','H mode @ fgr = 0.85', ...
             'Hybrid @ fgr = 0.85','Hybrid @ fgr = 0.39', ...
             'Steady State @ fgr = 0.86','Steady State @ fgr = 0.58'};
scenar = menu('JT60_SA scenario',scenario_liste);
switch scenar
case 1
	switch mode
	case 'cdr'
		Bt = 2.1;
		rb0 = sepa_option.ra .* Bt;
		ip  = 3e6; 

	otherwise
		Bt = 2.1;
		rb0 = sepa_option.ra .* Bt;
		ip  = 2.89e6; 

	end

	nbar  = 0.6 .* 1e20 .* (ip / 1e6) ./ (pi.* sepa_option.a .^ 2);
	plh   = 0;
	picrh = 0;
	pecrh = 0e6;
	pnbi  = 24e6;
	zeff  = 1.7;
	li    = 0.7;
	breakdown = 0.03;
	hmore = 1.0;
	xece  = 0.2;
	rtang = sepa_option.ra;
	einj  = 85e3;
	pnbi_p  = 24e6;
	pnbi_n  = 0e6;
	ane = 1;
	vane = 1.01;
	xiioxie = 1;
        fpped = 1;	   
        angle_nbi = 55;	
	l2hmul = 1;
	zext  =0;
	
case 2
	switch mode
	case 'cdr'
		Bt = 2.69;
		rb0 = sepa_option.ra .* Bt;
		ip  = 3.5e6; 

	otherwise
		Bt = 2.27;
		rb0 = sepa_option.ra .* Bt;
		ip  = 3.07e6; 

	end

	nbar  = 0.85 .* 1e20 .* (ip / 1e6) ./ (pi.* sepa_option.a .^ 2);
	plh   = 0;
	picrh = 0;
	pecrh = 0e6;
	pnbi  = 24e6;
	zeff  = 1.7;
	li    = 0.7;
	breakdown = 0.03;
	hmore = 1.0;
	xece  = 0.2;
	rtang = sepa_option.ra;
	einj  = 85e3;
	pnbi_p  = 24e6;
	pnbi_n  = 0e6;
	ane = 1;
	vane = 1.01;
	xiioxie = 1;
        fpped = 1;	   
        angle_nbi = 55;	
	l2hmul = 1;
	zext  =0;

case 3
	switch mode
	case 'cdr'
		Bt = 2.69;
		rb0 = sepa_option.ra .* Bt;
		ip  = 3.5e6; 

	otherwise
		Bt = 2.27;
		rb0 = sepa_option.ra .* Bt;
		ip  = 2.84e6; 

	end

	nbar  = 0.85 .* 1e20 .* (ip / 1e6) ./ (pi.* sepa_option.a .^ 2);
	plh   = 0;
	picrh = 0;
	pecrh = 4e6;
	pnbi  = 34e6;
	zeff  = 1.7;
	li    = 0.7;
	breakdown = 0.03;
	hmore = 1.11;
	xece  = 0.2;
	rtang = sepa_option.ra;
	einj  = 85e3;
	pnbi_p  = 24e6;
	pnbi_n  = 10e6;
	ane = 1;
	vane = 1.01;
	xiioxie = 1;
        fpped = 1;	   
        angle_nbi = 55;	
	l2hmul = 1;
	zext  =0.2;

case 4 
	switch mode
	case 'cdr'
		Bt = 2.68;
		rb0 = sepa_option.ra .* Bt;
		ip  = 5.5e6; 

	otherwise
		Bt = 2.26;
		rb0 = sepa_option.ra .* Bt;
		ip  = 5.5e6; 

	end

	nbar  = 0.39 .* 1e20 .* (ip / 1e6) ./ (pi.* sepa_option.a .^ 2);
	plh   = 0;
	picrh = 0;
	pecrh = 4e6;
	pnbi  = 34e6;
	zeff  = 1.7;
	li    = 0.7;
	breakdown = 0.03;
	hmore = 1.33;
	xece  = 0.6;
	rtang = sepa_option.ra;
	einj  = 85e3;
	pnbi_p  = 24e6;
	pnbi_n  = 10e6;
	ane = 1;
	vane = 1.01;
	xiioxie = 1;
        fpped = 1;	   
        angle_nbi = 55;	
	l2hmul = 1;
	zext  =0.3;

case 5 
	switch mode
	case 'cdr'
		Bt = 1.81;
		rb0 = sepa_option.ra .* Bt;
		ip  = 2.41e6; 

	otherwise
		Bt = 1.61;
		rb0 = sepa_option.ra .* Bt;
		ip  = 2.54e6; 

	end

	nbar  = 0.86 .* 1e20 .* (ip / 1e6) ./ (pi.* sepa_option.a .^ 2);
	plh   = 0;
	picrh = 0;
	pecrh = 7e6;
	pnbi  = 34e6;
	zeff  = 1.7;
	li    = 0.7;
	breakdown = 0.03;
	hmore = 1.33;
	xece  = 0.6;
	rtang = sepa_option.ra;
	einj  = 85e3;
	pnbi_p  = 24e6;
	pnbi_n  = 10e6;
	ane = 1;
	vane = 1.01;
	xiioxie = 1;
        fpped = 1;	   
	l2hmul = 1;
        angle_nbi = 55;	
	zext  =0.3;

case 6
	switch mode
	case 'cdr'
		Bt = 2.1;
		rb0 = sepa_option.ra .* Bt;
		ip  =2.75e6; 

	otherwise
		Bt = 1.86;
		rb0 = sepa_option.ra .* Bt;
		ip  = 2.87e6; 

	end

	nbar  = 0.58 .* 1e20 .* (ip / 1e6) ./ (pi.* sepa_option.a .^ 2);
	plh   = 0;
	picrh = 0;
	pecrh = 7e6;
	pnbi  = 34e6;
	zeff  = 1.7;
	li    = 0.7;
	breakdown = 0.03;
	hmore = 1.33;
	xece  = 0.6;
	rtang = sepa_option.ra;
	einj  = 85e3;
	pnbi_p  = 24e6;
	pnbi_n  = 10e6;
	ane = 1;
	vane = 1.01;
	xiioxie = 1;
        fpped = 1;	   
        angle_nbi = 55;	
	l2hmul = 1;
	zext  =0.3;

end

option_phy ={'Johner','Natural'};
op_phy = menu('Physical model',option_phy);

switch op_phy
case 1	
	op_ane = 1;
	op_scale = 0;
case 2
	op_ane = 0;
	op_scale = 6;
end

% point de depart du plasma
tt    = 0.1;
ipt   = 0.15e6;
flux = (16 - 3) .* 2 .* pi;
Rt    = 2.5;
at    = 0.5;
Kt    = 1;
dt    = 0;
lit   = 1.5;
vloop = 2 .* pi .* Rt .* 0.15;
flux  = flux - 2 * pi .* tt .* vloop;
breakdown = 0.03;
hmoret = 1;
xecet  = 0;
plht   = 0;
picrht = 0;
pecrht = 2e6;
pnbit  = 0e6;
pnbi_pt  = 0e6;
pnbi_nt  = 0e6;
xiioxie = 0;
dipdt = 0.17e6;
ngr     = 1e20 .* (ipt(end) /1e6) ./ (pi .* at(end) .^ 2);
nsat = min(ngr,0.06e20 .* (ipt(end) ./ 1e6) .* Rt(end) .* sqrt(gaz) ./ Kt(end) ./ at(end) .^ (5/2));
nbart = nsat ./ 2;

% plasma en fin de break down
tt    = cat(1,tt,4);
ipt   = cat(1,ipt,ipt + dipdt .* abs(diff(tt)));
flux = cat(1,flux,flux);
Rt    = cat(1,Rt,2.8);
at    = cat(1,at,0.7);
Kt    = cat(1,Kt,1);
dt    = cat(1,dt,0);
hmoret = cat(1,hmoret,1);
xecet  = cat(1,xecet,0);
plht   = cat(1,plht,0);
picrht = cat(1,picrht,0);
pecrht = cat(1,pecrht,7e6);
pnbit  = cat(1,pnbit,0);
pnbi_pt  = cat(1,pnbi_pt,0);
pnbi_nt  = cat(1,pnbi_nt,0);
dipdt = 0.4e6;
ngr     = 1e20 .* (ipt(end) /1e6) ./ (pi .* at(end) .^ 2);
nsat = min(nbar,min(ngr,0.06e20 .* (ipt(end) ./ 1e6) .* Rt(end) .* sqrt(gaz) ./ Kt(end) ./ at(end) .^ (5/2)));
nbart = cat(1,nbart,nsat ./ 2);

% plasma en debut de plateau
dip = ip - ipt(end);
dtt  = dip ./ dipdt;
tt  = cat(1,tt,tt(end) + dtt);
ipt = cat(1,ipt,ip);
flux = cat(1,flux,flux);
Rt    = cat(1,Rt,sepa_option.ra);
at    = cat(1,at,sepa_option.a);
Kt     = cat(1,Kt,min(sepa_option.zxup,sepa_option.zxdo));
dt     = cat(1,dt,0.8 .* min(sepa_option.rxup,sepa_option.rxdo));
hmoret = cat(1,hmoret,hmore);
xecet  = cat(1,xecet,xece);
plht   = cat(1,plht,plh);
picrht = cat(1,picrht,picrh);
pecrht = cat(1,pecrht,pecrh);
pnbit  = cat(1,pnbit,pnbi);
pnbi_pt  = cat(1,pnbi_pt,pnbi_p);
pnbi_nt  = cat(1,pnbi_nt,pnbi_n);
dipdt = 0;
sepa_option.ton = tt(end);
nbart = cat(1,nbart,nbar);

% fin plateau
tt  = cat(1,tt,tt(end) + 100);
ipt = cat(1,ipt,ip);
flux = cat(1,flux,flux);
Rt    = cat(1,Rt,sepa_option.ra);
at    = cat(1,at,sepa_option.a);
Kt     = cat(1,Kt,min(sepa_option.zxup,sepa_option.zxdo));
dt     = cat(1,dt,0.8 .* min(sepa_option.rxup,sepa_option.rxdo));
hmoret = cat(1,hmoret,hmore);
xecet  = cat(1,xecet,xece);
plht   = cat(1,plht,plh);
picrht = cat(1,picrht,0);
pecrht = cat(1,pecrht,pecrh);
pnbit  = cat(1,pnbit,0);
pnbi_pt  = cat(1,pnbi_pt,0);
pnbi_nt  = cat(1,pnbi_nt,0);
nbart = cat(1,nbart,nbar);

% arret chauffage
dipdt = 0.2e6;
tt  = cat(1,tt,tt(end) + 10);
ipt = cat(1,ipt,ipt(end) - (tt(end) -tt(end-1)) .* dipdt);
flux = cat(1,flux,flux);
Rt    = cat(1,Rt,sepa_option.ra);
at    = cat(1,at,sepa_option.a);
Kt     = cat(1,Kt,min(sepa_option.zxup,sepa_option.zxdo));
dt     = cat(1,dt,0.8 .* min(sepa_option.rxup,sepa_option.rxdo));
hmoret = cat(1,hmoret,hmore);
xecet  = cat(1,xecet,xece);
plht   = cat(1,plht,0);
picrht = cat(1,picrht,0);
pecrht = cat(1,pecrht,0);
pnbit  = cat(1,pnbit,0);
pnbi_pt  = cat(1,pnbi_pt,0);
pnbi_nt  = cat(1,pnbi_nt,0);
dipdt = 0;
sepa_option.toff = tt(end);
ngr     = 1e20 .* (ipt(end) /1e6) ./ (pi .* at(end) .^ 2);
nsat = min(nbar,min(ngr,0.06e20 .* (ipt(end) ./ 1e6) .* Rt(end) .* sqrt(gaz) ./ Kt(end) ./ at(end) .^ (5/2)));
nbart = cat(1,nbart,nsat);

% descente de courant
dipdt = 0.4e6;
dtt    = (ip(end) - 0.15e6) ./ dipdt;
tt  = cat(1,tt,tt(end) + dtt);
ipt = cat(1,ipt,0.15e6);
flux = cat(1,flux,flux);
Rt    = cat(1,Rt,2.5);
at    = cat(1,at,0.5);
Kt     = cat(1,Kt,1);
dt     = cat(1,dt,0);
hmoret = cat(1,hmoret,1);
xecet  = cat(1,xecet,xece);
plht   = cat(1,plht,0);
picrht = cat(1,picrht,0);
pecrht = cat(1,pecrht,0);
pnbit  = cat(1,pnbit,0);
pnbi_pt  = cat(1,pnbi_pt,0);
pnbi_nt  = cat(1,pnbi_nt,0);
dipdt = 0;
ngr     = 1e20 .* (ipt(end) /1e6) ./ (pi .* at(end) .^ 2);
nsat = min(ngr,0.06e20 .* (ipt(end) ./ 1e6) .* Rt(end) .* sqrt(gaz) ./ Kt(end) ./ at(end) .^ (5/2));
nbart = cat(1,nbart,nsat);

% vecteur temps
temps = linspace(min(tt),max(tt),301)';


z0dinput = zerod_scalaire(temps,rb0./Rt(1),Rt(1),at(1),Kt(1),dt(1),ip,nbar,plh,picrh,pecrh,pnbi,zeff,xece,hmore,0,0,flux(1), ...
                          'gaz',gaz,'ane',op_ane,'vane',vane,'l2hmul',l2hmul, ...
			  'scaling',op_scale,'xiioxie',xiioxie,'qdds',0.95,'kidds',3, ...
			  'zeff',0,'zmax',8,'zimp',6,'rimp',0.1,'matthews',1, ...
			  'angle_ece',90,'sens',1, ...
                          'angle_nbi',angle_nbi,'einj',einj,'rtang',rtang,'lhmode',3,'etalh',0.8, ...
			  'npar0',2,'freqlh',5,'wlh',0,'xlh',0.5,'dlh',0.5, ...
                          'fwcd',0, 'mino','H','cmin',0.1,'nphi',25,'freq',72, ...
                          'sitb',0,'tae',0,'smhd',0,'tmhd',inf,'rip',0,'li',li,'breakdown',breakdown);

z0dinput.cons.nbar = interp1(tt,nbart,temps,'linear');
z0dinput.cons.ip = interp1(tt,ipt,temps,'linear');
z0dinput.cons.pecrh = z0interp_cons(tt,pecrht,temps);
z0dinput.cons.xece = z0interp_cons(tt,xecet,temps);
z0dinput.cons.plh =  z0interp_cons(tt,plht,temps);
%  switch op_phy
%  case 1	
%  
%  	z0dinput.cons.hmore = interp1(tt,hmoret,temps,'linear');
%  	z0dinput.option.sitb = 0;
%  	z0dinput.option.smhd = 100;
%  	z0dinput.option.tmhd = 1000000;
%  case 2
%  	z0dinput.cons.hmore = ones(size(z0dinput.cons.hmore));
%  	z0dinput.option.sitb = 3;
%  	z0dinput.option.smhd = 0;
%  	z0dinput.option.tmhd = 10;
%  end
z0dinput.cons.hmore = interp1(tt,hmoret,temps,'linear');
z0dinput.geo.R = pchip(tt,Rt,temps);
z0dinput.geo.a = pchip(tt,at,temps);
z0dinput.geo.K = pchip(tt,Kt,temps);
z0dinput.geo.d = pchip(tt,dt,temps);
z0dinput.cons.zeff  = 1.7 + 2.3 .*  (min(z0dinput.cons.nbar) ./ z0dinput.cons.nbar) .^ 2.6;
indon = find((temps > sepa_option.ton) &(temps < sepa_option.toff));
z0dinput.cons.zeff(indon) = zeff;


% transition X point
z0dinput = z0separatrix(z0dinput,sepa_option);
z0dinput.geo.b0 = rb0 ./z0dinput.geo.R;

% chema chauffage
if pnbi_n == 0
	z0dinput.cons.pnbi = z0interp_cons(tt,pnbi_pt,temps);
	z0dinput.option.rtang = 0;
	z0dinput.option.angle_nbi = 0;

else
	z0dinput.cons.pnbi = z0interp_cons(tt,pnbi_nt,temps);
	z0dinput.cons.picrh = z0interp_cons(tt,pnbi_pt,temps);
	z0dinput.option.cmin = 0.05;
	z0dinput.option.freq = mean(15.2 .* z0dinput.geo.b0 .*  z0dinput.geo.R ./ (z0dinput.geo.R + 0.8 .* z0dinput.geo.a));
	z0dinput.option.einj = 5e5;
	z0dinput.option.zext = zext;
	z0dinput.option.rtang = 2.5;
end


[zs,infovoid,profli] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
%[zs,infovoid,profli] = zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);

post.z0dinput = z0dinput;
post.zerod    = zs;
post.profil0d =profli;
op0d =post.z0dinput.option;


indmes = find(z0dinput.cons.temps > 80,1);
rfan = (3.56e6 + 14.03e6) ./ 3.56e6 ;
padd = z0dinput.cons.picrh + z0dinput.cons.pecrh + z0dinput.cons.pnbi + zs.pohm + zs.plh;
betanth = zs.wth .* (1.6.*pi./3) .* z0dinput.geo.a ./ zs.vp ./ z0dinput.geo.b0 ./ zs.ip;
fprintf('b0 = %g, R = %g, a = %g, Ka = %g, da = %g, epsi = %g\n', z0dinput.geo.b0(end-2),  z0dinput.geo.R(end-2), ...
        z0dinput.geo.a(indmes),z0dinput.geo.K(indmes),z0dinput.geo.d(indmes),z0dinput.geo.R(indmes)./z0dinput.geo.a(indmes));
fprintf('Ip = %g, qa = %g, q95 = %g, q0 = %g, qmin = %g, li = %g\n', ...
        zs.ip(indmes)./1e6,zs.qa(indmes),zs.q95(indmes), ...
        zs.q0(indmes),zs.qmin(indmes),zs.li(indmes));
fprintf('nbar/ngr = %g,<nD> = %g, <nT> = %g, <nHe> =%g,<nimp> =%g\n',zs.nbar(indmes) ./1e20./  (zs.ip(indmes) / 1e6) .*  ...
        (pi.* z0dinput.geo.a(indmes) .^ 2),zs.nDm(indmes)./1e19,zs.nTm(indmes)./1e19,zs.nhem(indmes)./1e19,zs.nimpm(indmes))./1e19;
fprintf('taue = %g, tauhe = %g , HH = %g, Wth = %g, W = %g, Wfast = %g\n', ...
        zs.taue(indmes),zs.tauhe(indmes),zs.taue(indmes)./zs.tauh(indmes).*zs.hitb(indmes),  ...
	zs.wth(indmes)./1e6,zs.w(indmes)./1e6,(zs.w(indmes) -  zs.wth(indmes)) ./1e6);
fprintf('Plh = %g, Picrh = %g, Pecrh = %g, Pnbi = %g, Pline = %g, Pbrem = %g, Pcyclo = %g\n',zs.plh(indmes)/1e6, ...
        zs.picrh(indmes)/1e6,zs.pecrh(indmes)/1e6,zs.pnbi(indmes)/1e6,(zs.prad(indmes)+ zs.pradsol(indmes))/1e6,zs.pbrem(indmes)/1e6, zs.pcyclo(indmes) ./1e6);
fprintf('<ne> = %g, ne0 = %g, ne0/<ne> = %g, <Te> = %g , Te0 = %g, Te0/<Te> = %g,  Ti/Te = %g\n',zs.nem(indmes)/1e19,zs.ne0(indmes)/1e19, ...
         zs.ne0(indmes) ./ zs.nem(indmes) ,zs.tem(indmes)/1e3,zs.te0(indmes)/1e3,zs.te0(indmes) ./ zs.tem(indmes),zs.tite(end -2));
fprintf('Pfus_in = %g, Pfus = %g, Q = %g, betan(total) = %g , betan(th) = %g, frad = %g  \n',zs.pfus(indmes)/1e6,zs.pfus(indmes).* rfan./1e6, ...
         rfan .* zs.pfus(indmes) ./ padd(indmes),zs.betan(indmes).*100,betanth(indmes).*100, ...
	 (zs.prad(indmes) + zs.pbrem(indmes) + zs.pcyclo(indmes)) ./ zs.pin(indmes));
fprintf('Iboot/Ip = %g, Icd/Ip = %g, Vloop = %g , Zeff(line) = %g, li =%g\n',zs.iboot(indmes) ./zs.ipar(indmes),zs.icd(indmes) ./ zs.ipar(indmes), ...
        zs.vloop(indmes),zs.zeff(indmes)./zs.zmszl(indmes),zs.li(indmes));
fprintf('Vol = %g, Section = %g, Sext = %g , Llcms = %g\n',zs.vp(indmes),zs.sp(indmes),zs.sext(indmes),zs.peri(indmes));
fprintf('frHe = %g, frImp = %g @ Z= %d  &  %g @ Z = %d, <Zeff> = %g \n', ...
         zs.nhem(indmes)./zs.nem(indmes),zs.nimpm(indmes)./zs.nem(indmes),op0d.zimp, ...
	 zs.nimpm(indmes)./zs.nem(indmes).* op0d.rimp,op0d.zmax, ...
	 (zs.n1m(indmes) + 4 .* zs.nhem(indmes) +zs.nimpm(indmes) .* op0d.zimp .^ 2 + ...
	 zs.nimpm(indmes) .* op0d.zmax  .^ 2 .* op0d.rimp) ./ zs.nem(indmes) );


z0plotsc


% consigne positif/negatif
if pnbi_n == 0
	zs.ndd_nbi_thp = zs.ndd_nbi_th;
	zs.ndd_nbi_thn = 0.*z0dinput.cons.pnbi;
	pnbi_nt = 0.*z0dinput.cons.pnbi;
	pnbi_pt = z0dinput.cons.pnbi;
	% correction du beam-beam
	fprintf('BB pp:');
	zs.ndd_nbi_nbi_p  = 	z0beambeam(zs.temps,profli.xli,profli.vpr,profli.nep,profli.tep,profli.n1p,profli.tip,profli.zeff,profli.pnbi, ....
                            zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,zs.nem,85e3,85e3,pnbi_pt,pnbi_pt,zs.mu0_nbi,1);
			    
	fprintf('BB nn:');	
	zs.ndd_nbi_nbi_n   = 0 .* zs.ndd_nbi_nbi_p; 
	
	fprintf('BB pn:');
	zs.ndd_nbi_nbi_pn  = 0 .* zs.ndd_nbi_nbi_p;  
else
	pnbi_nt = z0dinput.cons.pnbi;
	pnbi_pt = z0dinput.cons.picrh;
	option_l = z0dinput.option;
	option_l.einj = 85e3;
	option_l.cmin = 0;
	[void1,void2,zs.ndd_nbi_thp] = z0neutron_dd(option_l,cons,zs,profli,pnbi_pt); 
	option_l.einj = 500e3;
	option_l.cmin = 0;
	[void1,void2,zs.ndd_nbi_thn] = z0neutron_dd(option_l,cons,zs,profli,pnbi_nt); 

	
	% correction du beam-beam
	fprintf('BB pp:');
	zs.ndd_nbi_nbi_p  = z0beambeam(zs.temps,profli.xli,profli.vpr,profli.nep,profli.tep,profli.n1p,profli.tip,profli.zeff,profli.picrh, ....
                            zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,zs.nem,85e3,85e3,pnbi_pt,pnbi_pt,0 .* zs.mu0_nbi,1);
	fprintf('BB nn:');
	zs.ndd_nbi_nbi_n  = z0beambeam(zs.temps,profli.xli,profli.vpr,profli.nep,profli.tep,profli.n1p,profli.tip,profli.zeff,profli.pnbi, ....
                            zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,zs.nem,500e3,500e3,pnbi_nt,pnbi_nt,zs.mu0_nbi,1);
	fprintf('BB pn:');
	zs.ndd_nbi_nbi_pn  = z0beambeam(zs.temps,profli.xli,profli.vpr,profli.nep,profli.tep,profli.n1p,profli.tip,profli.zeff,profli.pnbi, ....
                            zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,zs.nem,85e3,500e3,pnbi_pt,pnbi_nt,zs.mu0_nbi,0);

			    
end


%  fprintf('BB total:');
%  zs.ndd_nbi_nbi_tot  = z0beambeam_test(zs.temps,profli.xli,profli.vpr,profli.nep,profli.tep,profli.n1p,profli.tip,profli.zeff,profli.pnbi, ....
%                              zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,zs.nem,85e3,500e3,pnbi_p,pnbi_n,zs.mu0_nbi,1);

% update des donnees



zs.ndd_total = zs.ndd_th + zs.ndd_nbi_thp + zs.ndd_nbi_thn + zs.ndd_nbi_nbi_p + zs.ndd_nbi_nbi_n + zs.ndd_nbi_nbi_pn;

post.zerod    = zs;
post.z0dinput.config_jt60.shape = shape_liste{shape};
post.z0dinput.config_jt60.scenario = scenario_liste{scenar};
post.z0dinput.config_jt60.gaz      =  gaz_liste{gaz};
post.z0dinput.config_jt60.op_phy = option_phy{op_phy};


JT60_SA_plot;

