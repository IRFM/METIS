% script - d'appel de helena a partir du zerod pour les etudes de demo
% script - du  plot des energy 0d
% calcul - 0d
% script - a inclure dans zero1t
% fonction - pour l'initialisation du zerod en mode scalaire
% script - du  plot de l'equilibre 0d
% formule - de Sauter
% appel - de lhacc_sym pour plot
% attention - ici nip est la densite de HDT
% calcul - du piquage de la pression
% script - du  plot des temperature 0d
% interface - graphique pour z0dacces
% calcul   - 
%---WARNING--- la 1ere ligne de commentaire de ZBOOT0DIFF.M est vide
% constante - physique (phys)
% calcul - d'une aproximation de bpol
% test - de l'inflence du piquage de la source sur le piquage de pe
%  Z0DPROFINFO  - courte description  
% extraction - des donnees scalaires
% diagramme - LH
% utilisee - pour l'interraction faisceau plasma D(d,n)He3
% fraction - de l'energie deposee sur les ion
% cette fonction tabule <sigma*v> pour l'interraction faisceau (D) - plasma (D)
% calcul - de la fonction de distribution de NBI
% calcul - de la vitesse de rotation moyenne
% coefficient - de geometrie
% calcul - des moments
% tri - des structure et affichage dans l'ordre
% script - de plot de la temperature ionique
% Z0DMOVIE  - cree un film a partir des donnees de Metis
% cette - fonction calcule le volume du plasma sa section et la surface externe
% calcul - de la puissance de fusion
% calcul - le flux de neutron dd
% script - de test de l'electroneutarlite
% cette - fonction extrait des donnees cronos le donnees du 0d
% integartion rapide de dn/dt = - n/tau + s en mode discret
% disp('callback - : ')
% version - matricielle du calcul des section efficaces
% script - de run de demo
% script - du  plot de LH 0d
% script - du  plot des puissance 0d
% script - du  plot des temperature 0d
% section - efficace effective pour l'injection de neutre
% calcul - la puissance de fusion du a l'interaction faisceau-plasma (idn)
% section - efficace T(d,n)He4 toutes energies
% script - de run de demo
% fonction - assistant d'appel du 0d cronos
% diagramme - LH
% securite - sur zeff
% calcul - le contenu en helium du plamsa
% visulaisation - 3D des profils 0D
% dans - metis la dirction toroidal est dans le sens du courant plasma, 
% script - de test de zsauter0d vs nclass dans cronos
% cette - fonction donne la solution pour 1 appel
% script - du  plot des densites 0d
% Z0DXDT - calul la derivee temporelle (difference simple)
% script - du  plot des DDS pour TS
% le - fit se fait sur 80% du rayon
% script - de test du zerod
% calcul - rapide avec extraction de point pertinent
% on - travail a ip donne pour le 0D
% cette fonction tabule <sigma*v> pour l'interraction faisceau (T) - plasma (D)
% script - du  plot des puissance 0d
% calcul - du courant du au runaway pour le debut du plasma
%  ZCLEAR0DSEPA  - courte description  
% script - de run de JTC_SU
% diagramme - LH
% calcul - de l'bsorbtion de NBI pour le zerod
% calcul - 0d a partir des donnees de cronos
% script - du  plot des puissance 0d
% calcul - des profils de densite en mode H
% script - pour le test du calcul de li avec varoition de geometrie
% perte - ripple : fci = couche centrale, 10% H dans D par defaut
% script - de test de z0rot2
% visulaisation - 3D des profils 0D
%  NOM - DE LA FONCTION  courte description  
% calcul - du cout de la machine
% script - d'appel du zerod pour le fit de l'efficacite de LH
% fonction - de calcul de l'estimation des sources de neutrons DD 
% script - pour le plot des neutrons
% pour - test :
% calcul - le rayonnement 0d
% script - du  plot des confinement 0d
% calul - de l'integrant pour sigmavnbitplasmad
%---WARNING--- la 1ere ligne de commentaire de SCENARZERODFCI.M est vide
% fonction - bootstrap Sauter corrigee (j*b)
% test - 
% cette - fonction calcul les effets du chauffage minoriatire a l'harmonique 1 de ICRH
% script - de verification des lois d'echelles
% calcul - le flux de neutron dd
% script - pour le test de z0icrh
% parametres - physiques
% test - du calcul de wbp
% limite  - sur tau
% calcul - la puissance de fusion du a l'interaction faisceau-plasma (idn)
% cette - fonction extrait des donnees cronos le donnees du 0d
% nbar   - = nl /2/a (10^20 m^-3)
%  FSHEAR  - courte description  
% cette - fonction cree un jeu de donnees cronos a partir des donnees du  0d
% test - de gg0d
% calcul - un profil de q monotonic
% calcul - de l'effet du courant de retour
%  CLOSE_0PLOT  - courte description  
taus - = zpmean(temps,pinjo,taus);  % approximation : une seule constante de temps;
% script - du  plot de la geometrie du 0d
	rapfile = - sprintf('z0d_%s_%d',info.machine,fix(info.shot));
% input - npar0
% constantes - 
% fit - des mesure d'efficacite de generation de courant FWCD
% script - de comparaison dt taux de neutron pour TS
%  Z0LOGLIN  - courte description  
% script - du  plot de la convergence du 0d
% declaration - des parametres
% elagrissement - du depot du fait de la largeure d'orbite
% script - de run du rampup de ITER
% fonction - de calcul d'un profil a 3 parametres
% cette - fonction estime le champ electrique radial a partir des donnees du 0D
% script - de preparation d'un temps pour un simulation demo
% test - de l'inflence du piquage de la source sur le piquage de pe
% fonction - de morphing des surface de flux
% fonction - auxiliaire pour le calul de nhem
% test - de zsauter0d
