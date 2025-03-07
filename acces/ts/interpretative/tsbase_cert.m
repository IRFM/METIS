function [cert,vers,date,heure,unix,uniy,uniz,nomdon,type]=tsbase_cert(certdon)

% TSBASE_CERT	Decodage du champ certifications sur la donnée rendu par TSBASE.
%	[cert,vers,date,heure,unix,uniy,uniz,nomdon,type]=tsbase_cert(certdon)
%
% cette fonction decode le champ de renseignements certdon rendu par tsbase, 
% qui contient, entre autres, la certification de la donnee.
%
% cert	=  valeurs des certifications (une par signal).
% vers	=  numero de version du traitement qui a produit la donnee.
% date	=  date de l'execution du traitement (chaine de caractere).
% heure =  heure de l'execution du traitement (chaine de caractere).
% unix	=  unite de la coordonnee 1 (chaine de caractere).
% uniy	=  unite de la coordonnee 2 (chaine de caractere).
% uniz	=  unite de la coordonnee 3 (chaine de caractere).
% nomdon=  nom de la donnée lue (chaine de caractere)
% type  =  type de la donne lue ('D' ou 'T')
%
% NB: dans le vocabulaire de la base,la coordonnee 1 correspond a la mesure.
%    cf manuel d'utilisation de TSLib (note ts-xx-93-xx)
%  Modification  YB le 3/10/94 :
%  Prise en compte nom et type de la donnée lue YB & RM le 3/6/02
%    renseignement particulierement interessant pour une donnee generique
%  si la matrice d'entree est vide, on retourne des valeurs par defaut 
cert = -9999;
vers = -9999;
date = '01/01/1900';
heure = '00:00:00';
unix = ' ';
uniy = ' ';
uniz = ' ';
nomdon ='';
type ='';

if ~isempty(certdon)
  nbt =length(certdon);
  nb_cert =certdon(1);
  deb=2;fin=1+nb_cert;
  cert = certdon(deb:fin);
  deb=fin+1;fin=fin+1;
  vers =certdon(deb:fin);
  deb=fin+1;fin=fin+8;
  unix =char(certdon(deb:fin)');
  deb=fin+1;fin=fin+8;
  uniy = char(certdon(deb:fin)');
  if nbt > 2+nb_cert+16+10+8+21+1
     deb = fin+1;fin=fin+8;
     uniz = char(certdon(deb:fin)');
  end
  deb = fin+1;fin = fin+10;
  date = char(certdon(deb:fin)');
  deb = fin+1;fin=fin+8;
  heure = char(certdon(deb:fin)');
  deb = fin+1;fin=fin+21;
  nomdon=char(certdon(deb:fin)');
  deb = fin+1;fin=fin+1;
  type=char(certdon(deb:fin)');
end
