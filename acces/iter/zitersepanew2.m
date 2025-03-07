function  sepa = zitersepanew2(temps,option)

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
	valeur.b0        = 5.3;      % magnetic field at R0 
	valeur.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
	valeur.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
	
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
	type.nbp       = 'integer';     
	type.mode      = 'string';     
	
	borne.rxup      = [0.05,0.9];   
	borne.zxup      = [0.5 3];
	borne.apup      = [0,25];
	borne.amup      = [0,90];  
	borne.ra        = [5,7];   
	borne.za        = [0,2];   
	borne.a         = [0.1,2.2];     
	borne.rxdo      = [0.05,0.9];   
	borne.zxdo      = [0.5 3];
	borne.apdo      = [0,25];
	borne.amdo      = [0,90];  
	borne.b0        = [2,6];   
	borne.nbp     = [35,250];     
	borne.mode     = {'elliptical','ellipse + hyperbola','hyperbola + ellipse','2 semi ellpise'};     
	
	
	defaut.rxup      = 0.466;     
	defaut.zxup      = 1.687;    
	defaut.apup      = 0;    
	defaut.amup      = 0;   
	defaut.ra      = 6.2;      
	defaut.za      = 0.65;      
	defaut.a       = 2;        
	defaut.rxdo      = 0.568;     
	defaut.zxdo      = 2.001;    
	defaut.apdo      = 22.46;    
	defaut.amdo      = 67.92;   
	defaut.b0      = 5.3;      
	defaut.nbp     = 201;     
	defaut.mode     = 'elliptical';     
	
	info.rxup      = 'upper triangularity (minor radius unit), X point distance from magnetic axis';
	info.zxup      = 'upper altitude X point (minor radius unit)';
	info.apup      = 'upper separatrix angle (R,X)  (LFS, degrees)';
	info.amup      = 'upper separatrix angle (-R,X) (HFS, degrees)';
	info.ra      = 'major radius (m) [6.2]';
	info.za      = 'altitude of the plasma (m) [6.2]';
	info.a       = 'minor radius (m) [2]';
	info.rxdo      = 'lower triangularity (minor radius unit), X point distance from magnetic axis';
	info.zxdo      = 'lower altitude X point (minor radius unit)';
	info.apdo      = 'lower separatrix angle (R,X)  (LFS, degrees)';
	info.amdo      = 'lower separatrix angle (-R,X) (HFS, degrees)';
	info.b0      = 'toroidal magnetic field at 6.2 m (T)';
	info.nbp     = 'number of points for the separatrix (depends on equilibrium module) [201]';
	info.mode = 'nature of the magnetic surface'
	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	sortie.description = 'ITER separatix parameters';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	sepa = sortie;
	return
end



if nargin <2
	error('you need a time base and a data structure');
end
if isempty(temps)
	temps = 1;
end

% compatibilite
rxu  =	option.rxup;     
zxu  =	option.zxup;    
apu  =	option.apup;    
amu  =	option.amup;   
ra  =	option.ra ;      
a   =	option.a;        
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
end

if mode == 2
  lim2 = atan(0.5*geom(5)/(1-geom(6)))*180/pi;

  if geom(8) < lim2

    geom(8) = lim2+0.1;

  end
end
if mode == 3
  if geom(4) < lim4

    geom(4) = lim4+0.1;

  end
  lim8 = atan(0.5*geom(5)/(1-geom(6)))*180/pi;
  if geom(8) > lim8

    geom(8) = lim8-0.1;

  end
end

[x,y]=separatrice(geom,mode);
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
fprintf('separatrix volume @0.95 = %g m^3 \n',0.95 .* vv);


% correction du b0, calcul sur l'axe geometrique
b0 = 6.2 .* option.b0 ./ ra;

% calcul des parametres pour helena
hr0     = ra;
hz0     = option.za;
ha      = a;
he1     = geom(1);
htrl    = geom(2); 
htrh    = geom(6); 
za = option.za;
% plot si pas d'output
h = findobj(0,'type','figure','tag','sepa_iter');
if isempty(h)
       h=figure('tag','sepa_iter');
else
       figure(h);
end   
   
 clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times','color',[1 1 1])
 
digi = load(fullfile(fileparts(which('zitersepa')),'itershape.mat'));

plot(rr,zz+za,'bo',R,Z+za,'g',digi.R1,digi.Z1,':k',digi.R2,digi.Z2,'-.k');
axis([0 9 -4.5 4.5]);
axis('equal');
xlabel('R (m)')
ylabel('Z (m)')
title('ITER separatrix');

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
fprintf('Rmin = %g, Rmax = %g, Zmin = %g, Zmax = %g, \n Rzmin = %g, Rzmax = %g\n', ...
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

