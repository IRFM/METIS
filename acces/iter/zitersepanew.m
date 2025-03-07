function  separatrice = zitersepa(temps,option)

% declaration des parametres
if nargin <=1 
	valeur.rxup      = 0.4;     % decalage du point x haut p/r grand rayon(m)
	valeur.zxup      = 1.67;    % hauteur du point x haut (m)
	valeur.apup      = 0;    % angle de la separatrice au point x (R,X)  (cote LFS, degres)
	valeur.amup      = 0;    % angle de la separatrice au point x (-R,X) l'horizontale(cote HFS, degres)
	valeur.ra      = 6.2;      % grand rayon geometrique (m) [6.2]
	valeur.a       = 2;        % petit rayon geometrique (m) {2,1.85}
	valeur.rxdo      = 0.4;     % decalage du point x bas p/r grand rayon(m)
	valeur.zxdo      = 2;    % hauteur du point x bas (m)
	valeur.apdo      = 22.46;    % angle de la separatrice au point x bas (R,X)  (cote LFS, degres)
	valeur.amdo      = 67.92;    % angle de la separatrice au point x bas (-R,X) l'horizontale(cote HFS, degres)
	valeur.nbp     = 201;      % nombre de point decrivant la separatrice (selon module equilibre) [201]
	
	type.rx      = 'float';   
	type.zx      = 'float';
	type.ap      = 'float';
	type.am      = 'float';  
	type.ra      = 'float';   
	type.a       = 'float';     
	type.k       = 'float';  
	type.d       = 'float';   
	type.b0      = 'float';   
	type.nbp     = 'integer';     
	
	borne.rx      = [5.07,5.15];   
	borne.zx      = [-3.34,3.42];
	borne.ap      = [20,25];
	borne.am      = [65,70];  
	borne.ra      = [5,7];   
	borne.a       = [0.1,2.2];     
	borne.k       = [1,3];  
	borne.d       = [0,0.7];   
	borne.b0      = [0,5.3];   
	borne.nbp     = [35,250];     
	
	
	defaut.rx      = 5.09;     
	defaut.zx      = -3.36;    
	defaut.ap      = 22.46;    
	defaut.am      = 67.92;   
	defaut.ra      = 6.2;      
	defaut.a       = 2;        
	defaut.k       = 1.85;     
	defaut.d       = 0.49;     
	defaut.b0      = 5.3;      
	defaut.nbp     = 201;     
	
	info.rx      = 'rayon du point x (m)';
	info.zx      = 'hauteur du point x (m)';
	info.ap      = 'angle de la separatrice au point x /a l''horizontale (R,X)    (cote LFS, degres)';
	info.am      = 'angle de la separatrice au point x  /a l''horizontale (X,R)   (cote HFS, degres)';
	info.ra      = 'grand rayon geometrique (m) [6.2]';
	info.a       = 'petit rayon geometrique (m) {2,1.85}';
	info.k       = 'elongation (b/a) {1.85,1.97}';
	info.d       = 'triangularite moyenne (Cl + Cu) / 2*a {0.49,0.58}';
	info.b0      = 'champ toroidal a R = 6.2 m {5.3,5.18}';
	info.nbp     = 'nombre de point decrivant la separatrice (selon module equilibre) [201]';
	
	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	sortie.description = 'Parametres de la separatrice d''ITER';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	separatrice = sortie;
	return
end



if nargin <2
	error('il faut donnee une base-temps et un structure de parametres');
end
if isempty(temps)
	temps = 1;
end

% parametre de calcul des derivee :
l = 1e-2;

% compatibilite
rx  =	option.rx;     
zx  =	option.zx;    
ap  =	option.ap;    
am  =	option.am;   
ra  =	option.ra ;      
a   =	option.a;        
k   =	option.k;     
d   =	option.d;     
b0  =	option.b0;      
nbp = 	option.nbp ;     

% initialisation des variables de sortie
separatrice = [];


% calcul du z du point haut
zh  = min(4.2,zx + 2 .* k .* a);
za  = 0.5 .*  (zx + zh);
% max(Rmax) <= 8.2 m
% min(Rmin) >= 4.18 m
rmax = min(ra + a,8.2);
rmin = max(4.18,rmax  - 2 .* a);
ra   = 0.5 .* (rmin + rmax);
a    = 0.5 .* (rmax - rmin);

% r point haut
rh = 2 .* ra - rx  - 2 .* a .* d;

% point de contrainte sur la separatrice
rr = [];
zz = [];
% le point externe equatorial (Rmax)
rr(end+1) = rmax;
zz(end+1) = za - l ;
rr(end+1) = rmax;
zz(end+1) = za;
rr(end+1) = rmax - l /30;
zz(end+1) = za + l ;
% le point haut
rr(end+1) = rh + l;
zz(end+1) = zh - l/10;
rr(end+1) = rh;
zz(end+1) = zh;
rr(end+1) = rh - l;
zz(end+1) = zh - l/5;
indh       = length(zz); 
% le point interne (Rmin)
rr(end+1) = rmin ;
zz(end+1) = za + l ;
rr(end+1) = rmin ;
zz(end+1) = za;
rr(end+1) = rmin + l /30;
zz(end+1) = za - l ;
% le point avant le point x
rr(end+1) = rx - 1.7 .* l .* cos(am .* pi ./ 180);
zz(end+1) = zx + 2 .* l .* sin(am .* pi ./ 180);
rr(end+1) = rx - l  .* cos(am .* pi ./ 180);
zz(end+1) = zx + l  .* sin(am .* pi ./ 180);
% le point  x
rr(end+1) = rx;
zz(end+1) = zx;
% le point apres le point x
rr(end+1) = rx + l .* cos(ap .* pi ./ 180);
zz(end+1) = zx + l .* sin(ap .* pi ./ 180);
rr(end+1) = rx + 1.8 .* l .* cos(ap .* pi ./ 180);
zz(end+1) = zx + 2 .* l .* sin(ap .* pi ./ 180);
% le point externe equatorial (Rmax)
rr(end+1) = rmax - l /7;
zz(end+1) = za -l ;
rr(end+1) = rmax;
zz(end+1) = za ;
rr(end+1) = rmax;
zz(end+1) = za + l ;

% procedure d'optimisation de la forme
ok = 0;
zk = 2;
while (ok == 0)
	ok = 1;
	% angle
	alpha  = unwrap(angle((rr-ra)+ i .* (zz-za)));
	% reechantillonage par spline
	teta  = linspace(0,2*pi,nbp);
	try
	    R     = spline(alpha,rr,teta);
	    Z     = spline(alpha,zz,teta);
	catch    
	    hw    = warndlg('Il n'' a pas de solution avec ce jeux de parametres','Probleme lors du calcul de la separtrice'); 
	    waitfor(hw);
	    return
	end

	% rmin
	ind  = min(find((R < rmin)& (R == min(R))));
	if ~isempty(ind)
		ok =0;
		indi = min(find(alpha>teta(ind)));
		rrmem = rr;
		zzmem = zz;
		rr = cat(2,rrmem(1:(indi-1)),rmin + zk .* l,rrmem(indi:end));
		zz = cat(2,zzmem(1:(indi-1)),Z(ind),zzmem(indi:end));
                zk = zk /2;
	end
	ind  = min(find((R > rmax) & (R == max(R))));
	if ~isempty(ind) & (ok == 1)
		ok =0;
		indi = min(find(alpha>teta(ind)));
		rrmem = rr;
		zzmem = zz;
		rr = cat(2,rrmem(1:(indi-1)),rmax - l,rrmem(indi:end));
		zz = cat(2,zzmem(1:(indi-1)),Z(ind),zzmem(indi:end));	
	end		
	ind  = min(find((Z < zx) & (Z == min(Z))));
	if ~isempty(ind) & (ok == 1)
		ok =0;
		indi = min(find(alpha>teta(ind)));
		rrmem = rr;
		zzmem = zz;
		rr = cat(2,rrmem(1:(indi-1)),R(ind),rrmem(indi:end));
		zz = cat(2,zzmem(1:(indi-1)),zx + l,zzmem(indi:end));	
	end		
	ind  = min(find((Z > zh) & (Z == max(Z))));
	if ~isempty(ind) & (ok == 1)
		ok =0;
		indi = min(find(alpha>teta(ind)));
		rrmem = rr;
		zzmem = zz;
		rr = cat(2,rrmem(1:(indi-1)),R(ind),rrmem(indi:end));
		zz = cat(2,zzmem(1:(indi-1)),zh - l,zzmem(indi:end));	
	end		
end

% correction derivee point haut
%alpha  = unwrap(angle((rr-ra)+ i .* (zz-za)));
%rc     = abs((rr(indh)-ra)+ i .* (zz(indh)-za));
%alpha  = 2.* alpha(indh) - alpha(indh - 1);
%rr(indh) = ra + rc .* cos(alpha);
%zz(indh) = za + rc .* sin(alpha);

% angle
alpha  = unwrap(angle((rr-ra)+ i .* (zz-za)));
% reechantillonage par spline
teta  = linspace(0,2*pi,nbp);
R     = spline(alpha,rr,teta);
Z     = spline(alpha,zz,teta);

% calcul du volume
% abscisse R,  Z>za
indp = find(Z>=za);
rzp = R(indp);
zzp = Z(indp);
vp   = - 2 .* pi .* trapz(rzp,(zzp - za) .* rzp);
% abscisse R,  Z<za
indm = find(Z<=za);
rzm = R(indm);
zzm = Z(indm);
vm  = - 2 .* pi .* trapz(rzm,(zzm - za) .* rzm);
vv   = vp + vm;
fprintf('Volume separatrice @0.95 = %g m^3 \n',0.95 .* vv);


% correction du b0, calcul sur l'axe geometrique
b0 = 6.2 .* b0 ./ ra;

% calcul des parametres pour helena
hr0     = ra;
hz0     = za;
ha      = a;
he1     = k;
htrl    = asin((rx - ra) ./ a); 
htrh    = - asin((ra - 2 .* a.* d - rx) ./ a); 

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

plot(rr,zz,'bo',R,Z,'g',digi.R1,digi.Z1,':k',digi.R2,digi.Z2,'-.k');
axis([0 9 -4.5 4.5]);
axis('square');
xlabel('R (m)')
ylabel('Z (m)')
title('separatrice d''ITER');

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
separatrice.R    = R;       % vecteur R des points de la separatrice (m)
separatrice.Z    = Z;       % vecteur Z des points de la separtrice (m)
separatrice.r0   = hr0;     % grand rayon du plasma (m)
separatrice.z0   = hz0;     % centre geometricque du plasma en Z (m)
separatrice.a    = ha;      % petit rayon du plasma (m)
separatrice.e1   = he1;     % elongation (b/a)
separatrice.trl  = htrl;    % triangularite haute (definition entree de helena)
separatrice.trh  = htrh;    % triangularite basse (definition entree de helena)
separatrice.b0   = b0;      % champ toroidal a R = 6.2 m (T)
separatrice.vv   = vv;      % volume delimite par la separtrice (m^3)
separatrice.rr   = rr;      % vecteur R  des points definis avant spline (pour test)
separatrice.zz   = zz;      % vecteur Z des points definis avant spline  (pour test)

