% ZINITEQUI  valeurs par defaut du formulaire des donnees initiales de l'equilibre
%-----------------------------------------------------------------
% fichier : zinitequi.m 
% 
% fonction Matlab 5 : 
%	fonction qui definit les valeurs par defaut du formulaire "des donnees initiales de l'equilibre"
%  
% syntaxe :  
%   d_out = zinitequi
%  
% entrees :  
%  
% sorties :  
%   d_out  = structure contenant les valeurs par defaut
%  
% fonction �rite par J-F Artaud, poste 46-78
% version  1.7  du  29/09/2001  
%  
% liste des modifications :  
%  
%-------------------------------------------------------------------------------  

function sortie = zinitequi

if nargout > 0
	% declaration
	valeur.li          = 1;      % li de reference pendant la phase de chauffage [1]
	type.li            = 'float';
	borne.li           = [0.3,3];
	defaut.li          = 1;     
	info.li            = 'initial internal self [1]';
	valeur.ip          = 1.000000001e6;      % li de reference pendant la phase de chauffage [1]
	type.ip            = 'float';
	borne.ip           = [1e3,1e8];
	defaut.ip          = 1e3;
	info.ip            = 'initial plamsa current (A) [1]';

	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	sortie.description = 'initial data for equilibrium computation';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe

	return
end	

% lecture des donnees
geo   = evalin('base','data.geo');
geo.R = double(geo.R);
geo.Z = double(geo.Z);

ip   = evalin('base','data.cons.ip');
ipr  = ip; 
flux = evalin('base','data.cons.flux');
li   = evalin('base','data.exp.li');
t    = evalin('base','data.gene.temps');
kini = evalin('base','param.gene.kmin');
x    = evalin('base','param.gene.x');
mu0  = evalin('base','param.phys.mu0');


% dialogue pour fixe li :
info = zinitequi;
zassignin('base','initequi.li',li(kini));                  
zassignin('base','initequi.ip',ip(kini));                  
h    = zuicreefunform('zinitequi','initequi',1);
set(h,'name',sprintf('initial data for equilibrium computation'));
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end
%
li   = evalin('base','initequi.li');
ip   = evalin('base','initequi.ip');

%
if geo.mode(kini) == 2
	rmin  = min(geo.R(kini,:));
	rmax  = max(geo.R(kini,:));
	ra    = 0.5 .* (rmin + rmax);
	za    = geo.Z(kini,min(find(geo.R(kini,:) == rmax)));
	a     = 0.5 .* (rmax - rmin);
	zmin  = min(geo.Z(kini,:));
	zmax  = max(geo.Z(kini,:));
	b     = 0.5 .* (zmax - zmin);
	k     = b ./ a;
else
	ra    = geo.r0(kini);
	a     = geo.a(kini);
	k     = geo.e1(kini);

end

% piquage de j :
piqj = min(10,max(0.1,(exp(li)-1.65)./0.89));
% le profil de courant (initial)
jmoy  = (1 - x .^ 2 ) .^ piqj;
% normalisation approchee (sans importance pour la suite, renormlise par l'equilibre)
inorm = (2 .* pi .* a .^ 2 .* k) .* trapz(x,jmoy.*x);
j0    = ip ./ inorm;
jmoy  = j0 .* jmoy;

% psi initial
peri   = 4 .* a .* sqrt(k)  .* ellie(sqrt((k - 1) ./ k));
grhor  = (0.95 +(1.22 .* peri ./ (2 .* pi .* a)  -1) .* (x .^ 2)) ./ ra;
rhomax = peri ./2 ./ pi;
bpol0  = mu0 .* j0 .* a .* k ./ (1+ piqj) ./ 2 ;
bpol   = bpol0 .* ((1 - (1 - x .^ 2) .^ (piqj + 1)) ./ (x +eps));
bpol(1)= 0;
psi    = rhomax  .* cumtrapz(x,bpol ./ grhor);

% mise a la taille
vt     = ones(size(t));
ve     = ones(size(x));
psi    = vt * psi;
if all(isfinite(flux))
	psi    = psi(:,end) * ve - psi + flux * ve ; % verifier le sens ...
else
	psi    = psi(:,end) * ve - psi ; % verifier le sens ...
end
jmoy   = vt * jmoy;
% assignation
zassignin('base','data.prof.psi',psi);
zassignin('base','data.prof.jmoy',jmoy);
zassignin('base','data.equi.grhor',vt * grhor);
zassignin('base','data.equi.rhomax',vt .* rhomax);
zassignin('base','data.prof.bpolm',vt * bpol);


% recalcul du li
vv       =  4 * pi ^ 2 * ra * rhomax ^2 * k;
li        =  2 .* vv .* trapz(x,bpol .^2 .* x ) ./ ...
            ((mu0 .* ip ) .^ 2  .* ra);
zassignin('base','data.gene.li',vt .* li);
ipr(kini) =ip;
zassignin('base','data.cons.ip',ipr);


% calcul du ptot
pe   = evalin('base','data.prof.pe');
pion = evalin('base','data.prof.pion');
ptot = pe + pion;
zassignin('base','data.prof.ptot',ptot);
if any(~isfinite(ptot(kini,:)))
  warndlg('NaN in initial pressure profile -> Pe ou Pion uninitialized','init trouble');
end

zuisavenonok;
