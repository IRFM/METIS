% test de zimport
nom ='data.prof.zeff';
data.cons.idn = (NaN+i *NaN) .* data.cons.idn;
mode =1;
source = 1;
%filename ='/usr/drfc/cgc/matlab5/zineb/data/22530ts_courant_TiegalTe.mat.gz';
%filename ='tprof@28334.1';
 filename = 'test.ascii';
type = 1;
nom_temps = 'times';
%nom_espace = 'rhofit';
nom_espace = 2;
nom_data = 'tefit';
coordonnee = 2;
valeur_defaut = -pi;
positif =0;
zimport(nom,mode,source,filename,type,nom_temps,nom_espace,nom_data,coordonnee,valeur_defaut,positif)