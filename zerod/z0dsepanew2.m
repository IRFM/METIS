function  sepa = z0dsepanew2(temps,option,from_gui)

% declaration des parametres
if nargin <=1
	valeur.rxup      = 0.466;     % upper triangularity (minor radius unit)
	valeur.zxup      = 1.687;    % upper altitude X point (minor radius unit)
	valeur.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
	valeur.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
	valeur.ra        = 6.2;       % major radius R0 (m) [6.2]
	valeur.za        = 0.65;       % altitude of the magnetic axis (m) [0.9]
	valeur.a         = 2;         % minor radius (m) [2]
	valeur.rxdo      = 0.568;     % lower triangularity (minor radius unit)
	valeur.zxdo      = 2.001;       % lower altitude X point (minor radius unit)
	valeur.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
	valeur.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
	valeur.b0        = 13.6;      % magnetic field at R0
	valeur.delta     = 1;      %  
	valeur.update_b0 = 'off';      % 
	valeur.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
	valeur.mode      = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
        valeur.filename  = '';
	valeur.ton       = 0;
	valeur.toff      = Inf;
	
	type.rxup      = 'float';
	type.zxup      = 'float';
	type.apup      = 'float';
	type.amup      = 'float';  
	type.ra        = 'float';   
	type.za        = 'float';   
	type.a         = 'float';     
	type.rxdo      = 'float';   
	type.zxdo      = 'float';
	type.apdo      = 'float';
	type.amdo      = 'float';  
	type.b0        = 'float';
	type.delta     = 'float';
	type.update_b0 = 'string';
	type.nbp       = 'integer';
	type.mode      = 'string';     
	type.filename  = 'string';     
	type.ton       = 'float';
	type.toff      = 'float';
	
	borne.rxup      = [-2,2];
	borne.zxup      = [0.05,5];
	borne.apup      = [0,90];
	borne.amup      = [0,90];  
	borne.ra        = [0.1,100];
	borne.za        = [-10,10];
	borne.a         = [0.1,100];
	borne.rxdo      = [-2,2];
	borne.zxdo      = [0.05,5];
	borne.apdo      = [0,90];
	borne.amdo      = [0,90];  
	borne.b0        = [1e-4,100];
	borne.delta     = [0.1,100];
	borne.update_b0 = {'on','off'};
	borne.nbp       = [35,250];
	borne.mode      = {'elliptical','ellipse + hyperbola','hyperbola + ellipse','2 semi ellipse'};
        borne.filename  = '';
	borne.ton       = [0,Inf];
	borne.toff       = [0,Inf];

	
	defaut.rxup     = 0.466;     
	defaut.zxup     = 1.687;    
	defaut.apup     = 0;    
	defaut.amup     = 0;   
	defaut.ra       = 6.2;      
	defaut.za       = 0.65;      
	defaut.a        = 2;        
	defaut.rxdo     = 0.568;     
	defaut.zxdo     = 2.001;    
	defaut.apdo     = 22.46;    
	defaut.amdo     = 67.92;
	defaut.b0       = 13.6;
	defaut.delta    = 1;
	defaut.update_b0 = 'off';     
	defaut.nbp      = 201;
	defaut.mode     = 'elliptical';     
	defaut.filename = '';     
	defaut.ton      = 0;
	defaut.toff     = Inf;
	
	mode.rxup     = 'advanced';
	mode.zxup     = 'advanced';
	mode.apup     = 'advanced';
	mode.amup     = 'advanced';
	mode.ra       = 'advanced';
	mode.za       = 'advanced';
	mode.a        = 'advanced';
	mode.rxdo     = 'advanced';
	mode.zxdo     = 'advanced';
	mode.apdo     = 'advanced';
	mode.amdo     = 'advanced';
	mode.b0       = 'advanced';
	mode.delta    = 'advanced';
	mode.update_b0 = 'advanced';     
	mode.nbp      = 'advanced';
	mode.mode     = 'advanced';     
	
	info.rxup      = 'upper triangularity (minor radius unit), X point distance from magnetic axis';
	info.zxup      = 'upper altitude X point (minor radius unit)';
	info.apup      = 'upper separatrix angle (R,X)  (LFS, degrees)';
	info.amup      = 'upper separatrix angle (-R,X) (HFS, degrees)';
	info.ra        = 'major radius (m) [6.2]';
	info.za        = 'altitude of the plasma (m) [6.2]';
	info.a         = 'minor radius (m) [2]';
	info.rxdo      = 'lower triangularity (minor radius unit), X point distance from magnetic axis';
	info.zxdo      = 'lower altitude X point (minor radius unit)';
	info.apdo      = 'lower separatrix angle (R,X)  (LFS, degrees)';
	info.amdo      = 'lower separatrix angle (-R,X) (HFS, degrees)';
	info.b0        = 'maximal toroidal magnetic field on the conductor (T)';
	info.delta     = 'minimal distance between conductor an plasma (m)';
	info.update_b0 = 'if = on, update B0 accordingly to field on conductor, space between conductor a wall and R0 value';
	info.nbp       = 'number of points for the separatrix (depends on equilibrium module) [201]';
	info.mode      = 'nature of the magnetic surface';
	info.filename  = 'filename of matfile that contains the (R,Z) points of separatrix : \n variables must be vectors name R and Z \n if the filename field is not empty, uses the data from this matfile';
	info.ton       = 'first time for X point separatrix';
	info.toff      = 'last time for X point separatrix';
	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
    try
        if evalin('base','isfield(z0dinput,''sepa_option'')')
            %evalin('base','sepa_option=z0dinput.sepa_option');
        else
            %evalin('base','clear sepa_option');
            valeur.ton  = evalin('base','min(z0dinput.cons.temps)');
            valeur.toff = evalin('base','max(z0dinput.cons.temps)');
            zassignin('base','sepa_option',valeur);
        end
        tinter = cat(2,evalin('base','min(z0dinput.cons.temps)'),evalin('base','max(z0dinput.cons.temps)'));
    catch
            valeur.ton  = 0;
            valeur.toff = Inf;
            zassignin('base','sepa_option',valeur);
            tinter = cat(2,0,Inf);
    end
    
	borne.ton  = tinter;
	borne.toff = tinter;


	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.mode=mode;
	sortie.interface=interface;
	
	sortie.description = 'ITER separatix parameters';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	sepa = sortie;
	return
end

if isappdata(0,'GUI_GLOBAL_BASIC_MODE') && (getappdata(0,'GUI_GLOBAL_BASIC_MODE') == 1) && (nargin == 3)
    option.update_b0 = 'off';
end 
if nargin <2
	error('il faut donnee une base-temps et un structure de parametres');
end
if isempty(temps)
	temps = 1;
end


% compatibilite
rxu  =	option.rxup;     
zxu  =	option.zxup;
apu  =	option.apup;    
amu  =	option.amup;   
ra   =	option.ra ;      
a    =	option.a;        
rxd  =	option.rxdo;     
zxd  =	option.zxdo;    
apd  =	option.apdo;    
amd  =	option.amdo;   
nbp = 	option.nbp;     

% initialisation des variables de sortie
sepa = [];


geom(1) = zxu;
geom(2) = rxu;
geom(3) = apu;
geom(4) = amu;
geom(5) = zxd;
geom(6) = rxd;
geom(7) = apd;
geom(8) = amd;

if strcmp(option.mode,'elliptical')
  mode = 1;
end

if strcmp(option.mode,'ellipse + hyperbola')
  mode = 2;
end

if strcmp(option.mode,'hyperbola + ellipse')
  mode = 3;
end

if strcmp(option.mode,'2 semi ellpise')
  mode = 4;
end

if mode ==1 
  lim2 = atan(0.5*geom(5)/(1-geom(6)))*180/pi;
  if geom(8) > lim2

    geom(8) = lim2-0.1;

  end
  lim1 = atan(0.5*geom(1)/(1-geom(2)))*180/pi;
  if geom(4) > lim1

    geom(4) = lim1-0.1;

  end
end

if mode == 2
  lim2 = atan(0.5*geom(5)/(1-geom(6)))*180/pi;
  if geom(8) < lim2

    geom(8) = lim2+0.1;

  end
%    lim3 = atan(0.5*geom(1)/(1-geom(2))*180/pi;
%    if geom(4) < lim3
%  
%      geom(4) = lim3+0.1;
%  
%    end
end
if mode == 3
  lim4 = atan(0.5*geom(1)/(1-geom(2)))*180/pi;
  if geom(4) < lim4

    geom(4) = lim4+0.1;

  end
  lim8 = atan(0.5*geom(5)/(1-geom(6)))*180/pi;
  if geom(8) > lim8

    geom(8) = lim8-0.1;

  end
end
if isfield(option,'filename') && ~isempty(option.filename)
	try
	    sepadata = load(option.filename);
	catch
	    if nargin == 3
	        [filename, pathname] = uigetfile('*.mat', 'Select file for LCFS');
		if isequal(filename,0) || isequal(pathname,0)
			disp('LCFS creation cancelled')
			sepa = [];
			return
		else
			option.filename = fullfile(pathname,filename);
			sepadata = load(option.filename);			
		end
	    else
		disp(lasterr);	
		sepa = [];
		return
	    end
	end
	if isfield(sepadata,'R') && isfield(sepadata,'Z')
	    if (size(sepadata.R,1) > 1) && (size(sepadata.R,2) > 1)
		fprintf('R and Z variables must be vector and not matrix: LCFS creation cancelled\n');	
		sepa = [];
		return	    
	    end
	    r        = sepadata.R(:);
	    z        = sepadata.Z(:);
	    if length(r) ~= length(z)
		fprintf('R and Z variables have not the same length: LCFS creation cancelled\n');	
		sepa = [];
		return	    
	    end
	else 
	    fprintf('file %s not contain R and Z variables: LCFS creation cancelled\n',option.filename);	
	    sepa = [];
	    return
	end
	
%  	  % convex LCFS
%  	  KH = sort(unique(convhull(r,z)));
%  	  if (length(KH) ~= length(r))
%  	      index_full = 1:length(r);
%  	      Rsepa = r(KH);
%  	      Zsepa = z(KH);
%  	      r = interp1(KH,Rsepa,index_full,'linear');
%  	      z = interp1(KH,Zsepa,index_full,'linear');
%  	      indbad_lcfs = find(~isfinite(r) | ~isfinite(z));
%  	      if ~isempty(indbad_lcfs)
%  		  r(indbad_lcfs) = [];
%  		  z(indbad_lcfs) = [];
%  	      end
%  	  end
	
	
	%	 suppression des points doubles
	indbad = find((diff(r) == 0)  & (diff(z) == 0));
	nbw = 101;
	while ~isempty(indbad) && nbw > 0
	    r(indbad) = [];
	    z(indbad) = [];
	    nbw = nbw -1;
	    indbad = find((diff(r) == 0)  & (diff(z) == 0));
	end   
	% detection point x
	r_mem = r;
	z_mem = z;
	for l = 1:(length(r) - 5)
		[in,t,u,ri,zi,re,ze] = intersection_2_droites(r(l),r(l+2),z(l),z(l+2), ...
				r((l+3):(end-2)),r((l+5):end),z((l+3):(end-2)),z((l+5):end),eps);
		if any(in)
%  				figure(22);clf;plot([r(l),r(l+2)],[z(l),z(l+2)],'r', ...
%  				r((l+3):end),z((l+3):end),'b',ri(in~=0),zi(in~=0),'ok',re(in~=0),ze(in~=0),'xk');drawnow	

			indx  = find(in,1);
			if abs(t(indx) - 1) < 1e-9
				indmin = l + 2;
			elseif   abs(t(indx)) < 1e-9
				indmin = l;
			else
				indmin = l + 1;
			end
			if abs(u(indx) - 1) < 1e-9 
				indmax = indx + l + 4;
			elseif   abs(u(indx)) < 1e-9
				indmax = indx + l + 2;
			else
				indmax = indx + l + 3;
			end
			
			r=r(indmin:indmax);
			z=z(indmin:indmax);
			break;
		end
	end
        if length(r) < (length(r_mem) /3)
	  r = r_mem;
	  z = z_mem;
	  % echec separation des branches
	end
	r = r(:);
	z = z(:);

	r0   = (min(r) + max(r)) ./ 2;
	mask    = (r == max(r));
	z0   = sum(z .* mask) ./ max(1,sum(mask));
	cc   = (r - r0) + sqrt(-1) .* (z - z0);
	thc  = unwrap(angle(cc));
	thc(thc <0) = thc(thc<0) + 2 .* pi;
	rhoc = abs(cc);
	[thc,indc] = sort(thc);
	rhoc       = rhoc(indc);
	rhoc = cat(1,rhoc,rhoc,rhoc);
	thc = cat(1,thc -2.*pi,thc,thc+2.*pi);
	indnok = find(diff(thc)<=0);
	thc(indnok) =[];
	rhoc(indnok)  = [];
	teta  = linspace(0,2*pi,nbp);
	rho = spline(thc,rhoc,teta);
        R  = r0 + rho .* cos(teta);
        Z  = z0 + rho .* sin(teta);
	R(end) = (R(1)+ R(end)) ./ 2;
	R(1)   = R(end);
	Z(end) = (Z(1)+ Z(end)) ./ 2;
	Z(1)   = Z(end);
	% calcul du volume
	% abscisse R,  Z>za
	indp = find(Z>=0);
	rzp = R(indp);
	zzp = Z(indp);
	vp   = 2 .* pi .* trapz(rzp,(zzp) .* rzp);
	% abscisse R,  Z<za
	indm = find(Z<=0);
	rzm = R(indm);
	zzm = Z(indm);
	vm  = 2 .* pi .* trapz(rzm,(zzm) .* rzm);
	vv   = abs(vp) + abs(vm);
	fprintf('Volume separatrice @0.95 = %g m^3 \n',0.95 .* vv);
	
	% recalcul des parametres sur le vecteur final
	rmin  = min(R);
	rmax  = max(R);
	ra    = 0.5 .* (rmin + rmax);
	za    = Z(min(find(R == max(R))));
	a     = 0.5 .* (rmax - rmin);
	zmin  = min(Z);
	zmax  = max(Z);
	b     = 0.5 .* (zmax - zmin);
	k     = b ./ a;
	rzmax = R(min(find(Z == max(Z))));
	rzmin = R(min(find(Z == min(Z))));
	cl    = ra - rzmin;
	cu    = ra - rzmax;
	d     = (cl+cu) ./2 ./ a;
	rx    = rzmin;
	zx    = zmin;
	
	% correction du b0, calcul sur l'axe geometrique
	b0 = option.b0 .* (1 - a ./ ra - option.delta ./ ra);
	
	% calcul des parametres pour helena
	hr0     = ra;
	hz0     = za;
	ha      = a;
	he1     = b;
	htrl    = cl ./ a; 
	htrh    = cu ./ a; 

        rr = sepadata.R(:);
	zz = sepadata.Z(:);
elseif isappdata(0,'GUI_GLOBAL_BASIC_MODE') && (getappdata(0,'GUI_GLOBAL_BASIC_MODE') == 1) && (nargin == 3)
	sepa = [];
	return
else
	if isdeployed || ispc
	    [x,y]=separatrice_matfile(geom,mode);
	else
	    try
		  [x,y]=separatrice(geom,mode);
	    catch
		  [x,y]=separatrice_matfile(geom,mode);
	    end
	end
	rr = x*a+ra;
	zz = y*a;
	% angle
	alpha  = abs(unwrap(angle((rr-ra)+ i .* (zz))));
	% reechantillonage par spline
	teta  = linspace(0,2*pi,nbp);
	ind = find(diff(alpha) == 0);
	alpha(ind) = [];
	rr(ind) = [];
	zz(ind) = [];
	
	R     = spline(alpha',rr',teta);
	Z     = spline(alpha',zz',teta);
	
	R(end) = (R(1)+ R(end)) ./ 2;
	R(1)   = R(end);
	Z(end) = (Z(1)+ Z(end)) ./ 2;
	Z(1)   = Z(end);
	% calcul du volume
	% abscisse R,  Z>za
	indp = find(Z>=0);
	rzp = R(indp);
	zzp = Z(indp);
	vp   = 2 .* pi .* trapz(rzp,(zzp) .* rzp);
	% abscisse R,  Z<za
	indm = find(Z<=0);
	rzm = R(indm);
	zzm = Z(indm);
	vm  = 2 .* pi .* trapz(rzm,(zzm) .* rzm);
	vv   = vp + vm;
	fprintf('Volume separatrice @0.95 = %g m^3 \n',0.95 .* vv);
	
	
	% correction du b0, calcul sur l'axe geometrique
	b0 = option.b0 .* (1 - a ./ ra - option.delta ./ ra);
	
	% calcul des parametres pour helena
	hr0     = ra;
	hz0     = option.za;
	ha      = a;
	he1     = (geom(1) + geom(5)) / 2;
	htrl    = geom(2); 
	htrh    = geom(6); 
	za = option.za;
	Z  = Z + za;
end


% recalcul des parametres sur le vecteur final
rmin  = min(R);
rmax  = max(R);
ra    = 0.5 .* (rmin + rmax);
za    = Z(min(find(R == max(R))));
a     = 0.5 .* (rmax - rmin);
zmin  = min(Z);
zmax  = max(Z);
b     = 0.5 .* (zmax - zmin);
k     = b ./ a;
rzmax = R(min(find(Z == max(Z))));
rzmin = R(min(find(Z == min(Z))));
cl    = ra - rzmin;
cu    = ra - rzmax;
d     = (cl+cu) ./2 ./ a;
rx    = rzmin;
zx    = zmin;
fprintf('Rmin = %g, Rmax = %g, Zmin = %g, Zmax = %g, \nRzmin = %g, Rzmax = %g\n', ...
         rmin,rmax,zmin,zmax,rzmin,rzmax);
fprintf('Ra = %g, Za = %g, Rx = %g, Zx = %g\n' , ...
         ra,za,rx,zx);
fprintf('a = %g, K = %g, d = %g\n',a,k,d);
fprintf('b0 = %g\n',b0);


% mise a la dimension pour la base temps
if length(temps) > 1
	vt       = ones(length(temps),1);
	hr0      = vt .* hr0;
	hz0      = vt .* hz0;
	ha       = vt .* ha;
	he1      = vt .* he1;
	htrl     = vt .* htrl;
	htrh     = vt .* htrh;
	b0       = vt .* b0;
	vv       = vt .* vv;
	R        = vt * R(:)';
	Z        = vt * Z(:)';
end

%structure de sortie
sepa.R    = R;       % vecteur R des points de la separatrice (m)
sepa.Z    = Z;       % vecteur Z des points de la separtrice (m)
sepa.r0   = hr0;     % grand rayon du plasma (m)
sepa.z0   = hz0;     % centre geometricque du plasma en Z (m)
sepa.a    = ha;      % petit rayon du plasma (m)
sepa.e1   = he1;     % elongation (b/a)
sepa.trl  = htrl;    % triangularite haute (definition entree de helena)
sepa.trh  = htrh;    % triangularite basse (definition entree de helena)
sepa.b0   = b0;      % champ toroidal a R = 6.2 m (T)
sepa.vv   = vv;      % volume delimite par la separtrice (m^3)
sepa.rr   = rr;      % vecteur R  des points definis avant spline (pour test)
sepa.zz   = zz;      % vecteur Z des points definis avant spline  (pour test)
if isfield(option,'update_b0')
  sepa.update_b0    = option.update_b0;
else
  sepa.update_b0    = 'on';
end


% display information if from_gui
if nargin > 2
  option
end
