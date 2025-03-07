% script de preparation du rampup pour iter a partir d'un fichier existant
clear data param post
zuiload(fullfile(getappdata(0,'root'),'/acces/iter/rampup_model.mat.gz'));

% creation des donnees sur 80 s 
t_ref  = [1.5,	3.5,	7.33,	11.2,	15,	21.2,	27.4,	33.6,	39.8,	46,	52.1,	58.3,	64.5,	70.7,	80];
ip_ref = [0.5,	1.5,	2.5,	3.5,	4.5,	5.5,	6.5,	7.5,	8.5,	9.5,	10.5,	11.5,	12.5,	13.5,	15];

% base temps
tmax = 80;
data.gene.temps = linspace(1.5,tmax,length(data.gene.temps))';

% consigne de ip 
data.cons.ip = interp1(t_ref,ip_ref,data.gene.temps,'linear','extrap') .* 1e6;

% separtrice
[Rsepa,Zsepa] = separatrice_scenario_2(data.cons.ip);
data.geo.R    = Rsepa;
data.geo.Z    = Zsepa;
data.geo.mode(:) = 2;
rb0           = data.geo.b0 .* data.geo.r0;
% calcul des moments
% recalcul des parametres sur le vecteur final
rmin  = min(Rsepa,[],2);
rmax  = max(Rsepa,[],2);
data.geo.r0    = 0.5 .* (rmin + rmax);
data.geo.a     = 0.5 .* (rmax - rmin);
mask           = (Rsepa == (rmax * ones(1,size(Rsepa,2))));
data.geo.z0    = sum(Zsepa .* mask,2) ./ max(1,sum(mask,2));
zmin  = min(Zsepa,[],2);
zmax  = max(Zsepa,[],2);
data.geo.e1    = 0.5 .* (zmax - zmin) ./ data.geo.a;
mask  = (Zsepa == (zmax * ones(1,size(Zsepa,2))));
rzmax = sum(Rsepa .* mask,2) ./ max(1,sum(mask,2));
mask  = (Zsepa == (zmin * ones(1,size(Zsepa,2))));
rzmin = sum(Rsepa .* mask,2) ./ max(1,sum(mask,2));
data.geo.trh1  = (data.geo.r0 - rzmax) ./ data.geo.a ;
data.geo.trb1  = (data.geo.r0 - rzmin) ./ data.geo.a ;
data.geo.ind = zeros(length(data.gene.temps),1);
data.geo.b0    = rb0 ./ data.geo.r0;

% densite
ne0 = (data.gene.temps - 1.5) .* (0.4 - 0.025) ./ (80 -1.5) + 0.025;
ne0 = ne0 .* 1e20;
fedge = 0.4;
data.prof.ne = ne0 * ((1 - fedge) * (1 - param.gene.x .^ 3) + fedge); 
% zeff 
tz = [1.5,10,20,30,40,50,80];
zz = [4,2.5,2,1.8,1.75,1.73,1.70];
data.cons.zeffm  = interp1(tz,zz,data.gene.temps,'linear','extrap');
% module d'impurete
% attention on n'a pas le transport de He
% attention au parametre zinebcompo
param.cons.impur.zeff = 1;
param.cons.impur.exposant= 0
param.cons.impur.cmin1= 1
param.cons.impur.cmin2= 0.0500
param.cons.impur.rimp= 0.0600
param.cons.impur.neoimp= 0
param.cons.impur.nhndonoff= 0
param.cons.impur.norme= 1
param.cons.impur.crad =1
param.cons.impur.lrad =0.1000
param.cons.impur.frad =0.7000

% calcul de ni 
[ae,nion,nmin1,nmin2,nimp1,nimp2]=zitercompoinit(data.prof.ne,data.cons.zeffm * ones(size(param.gene.x)), ...
						 param.cons.impur.cmin1,param.cons.impur.cmin2,param.cons.impur.rimp, ...
                                                 param.compo.z(1),param.compo.z(2),param.compo.z(3),param.compo.z(4),param.compo.z(5));
data.prof.ae = ae;
data.prof.ni = data.prof.ne .* ae;
data.impur.impur = cat(3,nion,nmin1,nmin2,nimp1,nimp2);
	
					
% temperature de bord
data.cons.te1 = 25 .* ones(size(data.gene.temps));
data.cons.ti1 = 25 .* ones(size(data.gene.temps));

% calcul de la pression initiale
data.prof.te(1,:) =  (500 -  data.cons.te1(1)) .* (data.prof.te(1,:) - data.prof.te(1,end)) ./  ...
                     (data.prof.te(1,1) - data.prof.te(1,end)) +  data.cons.te1(1);
data.prof.ti(1,:) =  (250 -  data.cons.ti1(1)) .* (data.prof.ti(1,:) - data.prof.ti(1,end)) ./  ...
                     (data.prof.ti(1,1) - data.prof.ti(1,end)) +  data.cons.ti1(1);

data.prof.te   =	 ones(size(data.gene.temps)) * data.prof.te(1,:);	     
data.prof.pe   =	 data.prof.ne .* data.prof.te .* param.phys.e; 	      
data.prof.ti   =	 ones(size(data.gene.temps)) * data.prof.ti(1,:);	     
data.prof.pion =	 data.prof.ni .* data.prof.ti .* param.phys.e; 	      


% configuration interne		     
param.gene.cn = 0.3;
param.gene.nbforce =1;



% fixe le li de depart
zinitequi