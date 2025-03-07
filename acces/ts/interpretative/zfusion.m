% ZFUSION calcul de la puissance et des produits de fusion 
%-------------------------------------------------------------------------
% fichier zfusion.m ->  zfusion
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule la puissance de chauffage du aux alpha et
% les sources de matiere du aux produits de fusion. 
%  
% syntaxe  :
%  
%     [source,matiere,memoire] = zfusion(cons,source,matiere,geo,equi,gazbord,prof,neo, ...
%                                impur,phys,compo,gene,memoire);
%
% entrees :
%
%     cons    =  struture des consignes de calcul de la fonction
%     source  =  modele de source generique, valeurs a 0 
%     matiere =  datak.source.n.fus (source d'impuretes du aux reactions de fusion) 
%     geo     =  datak.geo  (geometrie )
%     equi    =  datak.equi  (equilibre)
%     gazbord =  datak.cons.c (pas utiliser dans cette fonction, reserve pour d'autres modules)
%     prof    =  datak.prof  (profils)
%     neo     =  datak.neo  (pas utiliser dans cette fonction, reserve pour d'autres modules)
%     impur   =  datak.impur (zeff , profil de densite de gaz et d'impuretes ...)
%     phys    =  param.phys
%     compo   =  param.compo (numero atomique et nombre de masses des gaz)
%     gene    =  param.gene
%     memoire =  datak.memoire.fus (valeur de reference pour le dernier calcul complet, 
%                pas utiliser dans cette fonction, reserve pour d'autres modules)
%
% sorties :
% 
%     source  =  structure source completees 
%     matiere =  tableau des sources de matiere du aux reaction de fusion complete
%     memoire =  datak.memoire.fus (valeur de reference pour le dernier calcul complet, 
%                         pas utiliser dans cette fonction, reserve pour d'autres modules)
% 
% parametres : aucun
% 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.2, du 25/11/2003.
% 
% 
% liste des modifications : 
%
%   * 02/10/2001 -> ajout de la structure memoire en sortie
%   * 06/11/2001 -> ajout des autres sources de matieres et de pussances
%   * 18/11/2002 -> gros bug : manque facteur 1/2 dans taux de reactions
%   * 25/11/2003 -> gros bug de copier coller : un facteur 2 manque dans la 
%                   taux de recation D-T
%
%--------------------------------------------------------------
%
function [source,matiere,memoire]=zfusion(cons,source,matiere,geo,equi,gazbord,prof,neo,impur, ...
                                  phys,compo,gene,memoire)
                                  
% mode initialisation 
% fonction auto declarante                             
if nargin <=1 

	valeur     = [];            % pas de parametres
	type       = [];
	borne      = [];
	defaut     = [];
	info       = [];
	
	interface.ts        = '';   % pas d'interface  a definir (nom de la fonction d'interfacage avec les donnees TS)
	interface.jet       = '';   % pas d'interface  a definir (nom de la fonction d'interfacage avec les donnees Jet)
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;

	sortie.description = 'Calcul simple des produits de fuison et du chauffage du a alpha';   % description (une ligne) de la fonction
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	source=sortie;
	
	return

end

% debut de calcul
%initialisation de martiere 
matiere =zeros(size(matiere));

% informations sur les reactions                                  
[dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion(prof.ti);
                                  
% recherche de la presence des gaz
ind_p    = find(compo.z == 1 & compo.a == 1);
ind_d    = find(compo.z == 1 & compo.a == 2);
ind_t    = find(compo.z == 1 & compo.a == 3);
ind_he3  = find(compo.z == 2 & compo.a == 3);
ind_he4  = find(compo.z == 2 & compo.a == 4);

% % mise  zeros des source de matiere definie
% if ~isempty(ind_p)
% 	matiere(1,:,ind_p)= 0 .* matiere(1,:,ind_p);
% end
% if ~isempty(ind_d)
% 	matiere(1,:,ind_d)= 0 .* matiere(1,:,ind_d);
% end
% if ~isempty(ind_t)
% 	matiere(1,:,ind_t)= 0 .* matiere(1,:,ind_t);
% end
% if ~isempty(ind_he3)
% 	matiere(1,:,ind_he3)= 0 .* matiere(1,:,ind_he3);
% end
% if ~isempty(ind_he4)
% 	matiere(1,:,ind_he4)= 0 .* matiere(1,:,ind_he4);
% end

% fraction de la puissance deposee sur les ions
% reference : D. Sigmar and J. Joyce, Nuclear fusion, 11(3), 1971.
% formule etablie pour DT, (fit a partir d'une figure)
%  0.2 < te < 80 keV 
t=[0,10,20,30,40,50,60,70,80,100].*1e3;
f=[0.01,0.146,0.269,0.377,0.473,0.55,0.615,0.665,0.698,0.72];
fion = interp1(t,f,prof.te,'spline');

% le nombre de reaction par unite de volume et par seconde  est 
%    -> n1*n2*sigmav pour deux especes differentes (DT) et
%    -> n^2/2 *sigmav pour une espece (DD)
% cf : La fusion thermonucleaire, p 17

% la fonction traite couple par couple 
% 1- D + T  
if ~isempty(ind_d) & ~isempty(ind_t)
	% nombre de reaction par m^3 et par s
        % attention bug dans la formule :  especes differentes
	%nr_dt = 0.5 .* impur.impur(1,:,ind_d) .* impur.impur(1,:,ind_t) .* dt.sv;
	nr_dt = impur.impur(1,:,ind_d) .* impur.impur(1,:,ind_t) .* dt.sv;
	% le flux de neutrons
	source.neutron.dt =nr_dt;
	% source d'he4 du a DT
	if ~isempty(ind_he4)
		matiere(1,:,ind_he4)= matiere(1,:,ind_he4) + nr_dt;
	end
	% disparistion de D et T
	matiere(1,:,ind_d)= matiere(1,:,ind_d) - nr_dt;
	matiere(1,:,ind_t)= matiere(1,:,ind_t) - nr_dt;
	% densite de puissance deposee sur les ions
	source.ion = dt.he4 .* phys.e .* nr_dt .* fion;
	% densite de puissance deposee sur les electrons
	source.el = dt.he4 .* phys.e .* nr_dt .* (1 - fion);
end

% 2- D + D canal p (on utilise la meme repartition que pour dt)
if ~isempty(ind_d) 
	% nombre de reaction par m^3 et par s
	nr_ddp = 0.5 .* impur.impur(1,:,ind_d) .^2 .* dd_p.sv;
	% le flux de protons
	source.proton.dd =nr_ddp;
	% source de T du a DD
	if ~isempty(ind_t)
		matiere(1,:,ind_t)= matiere(1,:,ind_t) + nr_ddp;
	end
	% source de p du a DD
	if ~isempty(ind_p)
		matiere(1,:,ind_p)= matiere(1,:,ind_p) + nr_ddp;
	end
	% disparistion de D 
	matiere(1,:,ind_d)= matiere(1,:,ind_d) - nr_ddp .* 2;
	% densite de puissance deposee sur les ions
	source.ion = source.ion + dd_p.p .* phys.e .* nr_ddp .* fion;
	source.ion = source.ion + dd_p.t .* phys.e .* nr_ddp .* fion;
	% densite de puissance deposee sur les electrons
	source.el  = source.el  + dd_p.p .* phys.e .* nr_ddp .* (1 - fion);
	source.el  = source.el  + dd_p.t .* phys.e .* nr_ddp .* (1 - fion);
	
end

% 3- D + D canal n (on utilise la meme repartition que pour dt)
if ~isempty(ind_d) 
	% nombre de reaction par m^3 et par s
	nr_ddn = 0.5 .* impur.impur(1,:,ind_d) .^2 .* dd_n.sv;
	% le flux de neutron
	source.neutron.dd =nr_ddn;
	% source de He3 du a DD
	if ~isempty(ind_he3)
		matiere(1,:,ind_he3)= matiere(1,:,ind_he3) + nr_ddn;
	end
	% disparistion de D 
	matiere(1,:,ind_d)= matiere(1,:,ind_d) - nr_ddn .* 2;
	% densite de puissance deposee sur les ions
	source.ion = source.ion + dd_n.he3 .* phys.e .* nr_ddn .* fion;
	% densite de puissance deposee sur les electrons
	source.el  = source.el  + dd_n.he3 .* phys.e .* nr_ddn .* (1 - fion);
end

% 4- T + T  (on neglige la puissance porduite par cette reaction)
if ~isempty(ind_t) 
	% nombre de reaction par m^3 et par s
	nr_tt = 0.5 .* impur.impur(1,:,ind_t) .^2 .* tt.sv;
	% le flux de neutron
	source.neutron.tt =2 .* nr_tt;
	% source de He4 du a TT
	if ~isempty(ind_he4)
		matiere(1,:,ind_he4)= matiere(1,:,ind_he4) + nr_tt;
	end
	% disparition de T 
	matiere(1,:,ind_t)= matiere(1,:,ind_t) - nr_tt .* 2;
end

% 5- D + He3  (on utilise la meme repartition que pour dt)
if ~isempty(ind_d) & ~isempty(ind_he3)
	% nombre de reaction par m^3 et par s
        % attention bug dans la formule :  especes differentes
	%nr_dhe3 = 0.5 .* impur.impur(1,:,ind_d) .* impur.impur(1,:,ind_he3) .* dhe3.sv;
	nr_dhe3 = impur.impur(1,:,ind_d) .* impur.impur(1,:,ind_he3) .* dhe3.sv;
	% le flux de protons
	source.proton.dhe3 =nr_dhe3;
	% source d'He4 du a D + HE3
	if ~isempty(ind_he4)
		matiere(1,:,ind_he4) = matiere(1,:,ind_he4) + nr_dhe3;
	end
	% source de H du a D + HE3
	if ~isempty(ind_p)
		matiere(1,:,ind_p) = matiere(1,:,ind_p) + nr_dhe3;
	end
	% disparistion de D et He3
	matiere(1,:,ind_d)= matiere(1,:,ind_d) - nr_dhe3 ;
	matiere(1,:,ind_he3)= matiere(1,:,ind_he3) - nr_dhe3 ;
	% densite de puissance deposee sur les ions
	source.ion = source.ion + dhe3.p .* phys.e .* nr_dhe3 .* fion;
	source.ion = source.ion + dhe3.he4 .* phys.e .* nr_dhe3 .* fion;
	% densite de puissance deposee sur les electrons
	source.el  = source.el  + dhe3.p .* phys.e .* nr_dhe3 .* (1 - fion);
	source.el  = source.el  + dhe3.he4 .* phys.e .* nr_dhe3 .* (1 - fion);
end

% 6- T + He3  canal p + n (on utilise la meme repartition que pour dt)
if ~isempty(ind_t) & ~isempty(ind_he3)
	% nombre de reaction par m^3 et par s
        % attention bug dans la formule :  especes differentes
	%nr_the3pn = 0.5 .* impur.impur(1,:,ind_t) .* impur.impur(1,:,ind_he3) .* the3_pn.sv;
	nr_the3pn = impur.impur(1,:,ind_t) .* impur.impur(1,:,ind_he3) .* the3_pn.sv;
	% le flux de protons
	source.proton.the3  = nr_the3pn;
	% le flux de neutrons
	source.neutron.the3 = nr_the3pn;
	% source d'He4 du a T + HE3
	if ~isempty(ind_he4)
		matiere(1,:,ind_he4) = matiere(1,:,ind_he4) + nr_the3pn;
        end
	% source de H du a T + HE3
	if ~isempty(ind_p)
		matiere(1,:,ind_p) = matiere(1,:,ind_p) + nr_the3pn;
        end
	% disparistion de T et He3
	matiere(1,:,ind_t)= matiere(1,:,ind_t) - nr_the3pn ;
	matiere(1,:,ind_he3)= matiere(1,:,ind_he3) - nr_the3pn ;
	% densite de puissance deposee sur les ions
	source.ion = source.ion + the3_pn.p .* phys.e .* nr_the3pn .* fion;
	source.ion = source.ion + the3_pn.he4 .* phys.e .* nr_the3pn .* fion;
	% densite de puissance deposee sur les electrons
	source.el  = source.el  + the3_pn.p .* phys.e .* nr_the3pn .* (1 - fion);
	source.el  = source.el  + the3_pn.he4 .* phys.e .* nr_the3pn .* (1 - fion);
end

% 7- T + He3  canal D (on utilise la meme repartition que pour dt)
if ~isempty(ind_t) & ~isempty(ind_he3)
	% nombre de reaction par m^3 et par s
        % attention bug dans la formule :  especes differentes
	%nr_the3d = 0.5 .* impur.impur(1,:,ind_t) .* impur.impur(1,:,ind_he3) .* the3_d.sv;
	nr_the3d = impur.impur(1,:,ind_t) .* impur.impur(1,:,ind_he3) .* the3_d.sv;
	% source d'He4 du a T + HE3
	if ~isempty(ind_he4)
		matiere(1,:,ind_he4) = matiere(1,:,ind_he4) + nr_the3d;
        end
	% source de D du a T + HE3
	if ~isempty(ind_d)
		matiere(1,:,ind_d) = matiere(1,:,ind_d) + nr_the3d;
        end
	% disparistion de T et He3
	matiere(1,:,ind_t)   = matiere(1,:,ind_t) - nr_the3d ;
	matiere(1,:,ind_he3) = matiere(1,:,ind_he3) - nr_the3d ;
	% densite de puissance deposee sur les ions
	source.ion = source.ion + the3_d.p .* phys.e .* nr_the3d .* fion;
	source.ion = source.ion + the3_d.he4 .* phys.e .* nr_the3d .* fion;
	% densite de puissance deposee sur les electrons
	source.el  = source.el  + the3_d.p .* phys.e .* nr_the3d .* (1 - fion);
	source.el  = source.el  + the3_d.he4 .* phys.e .* nr_the3d .* (1 - fion);
end

