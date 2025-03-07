% GCRONOS traitement automatique cronos (diffusion du courant)
% ------------------------------------------------------------
% fonction Matlab5 : gcronos.m -> gcronos
%
% Cette fonction lance le cronos en mode automatique apres choc.
% Les donnees de tprof doivent etre disponibles.
%
% syntaxe :
%
%  description :
%
%    [equi,courant,transport,generic,polar,ece,mse,resume,scronpar,com] = ...
%                            gcronos(numchoc,{occurrence,option_str,fichier})
%
%  test de lecture des donnees dans la base :
% 
%       [equi,courant,transport,generic,polar,ece,mse,resume,scronpar,com] = ...
%                               gcronos(numchoc,occurrence,'noexec');
%
%
%  lecture de donnees dans un fichier en interactif : 
%
%       [equi,courant,transport,generic,polar,ece,mse,resume,scronpar,com] = ...
%                               gcronos(numchoc,occurrence,'',-1);
%
%  lecture de donnees dans un fichier (nom du fichier donne) : 
%
%       [equi,courant,transport,generic,polar,ece,mse,resume,scronpar,com] = ...
%                               gcronos(numchoc,occurrence,'','nom_du_fichier_avec_le_chemin','nom_au_format_batch');
%
%  lecture des donnees dans la base + execution :
% 
%       [equi,courant,transport,generic,polar,ece,mse,resume,scronpar,com] = ...
%                               gcronos(numchoc,occurrence);
%
%  utilisation du traitement comme commande batch sur hercule :
%
%        gcronos(numchoc,occurence,'batch');
%
%
% entrees :
%    numchoc    = numero du choc (entier)
%    occurrence = numero d'occurrence (entier, si possible le meme que tprof, defaut = 0)
%    option_str = chaine de carateres des options envoyee par le sequenceur (reserve)
%    fichier    = nom du fichier a charger. Ci ce champ est rempli, cronos n'est pas execute.
%                 le fichier est charge et les donnees sont retournees
% 
% sorties :
%
%    equi        =  strutcures des donnees relatives a l'equilibre a ecrire dans la base
%    courant     =  strutcures des donnees relatives a la diffusion du a ecrire dans la base
%    transport   =  strutcures des donnees relatives au transport de la chaleur a ecrire dans la base
%    generic     =  strutcures des donnees complementaires, utiles a l'analyse,a ecrire dans la base
%    polar       =  strutcures des donnees recalculer de la polarimetrie, a ecrire dans la base
%    ece         =  strutcures des donnees recalculer de l'ece , a ecrire dans la base
%    mse         =  strutcures des donnees recalculer de la mse , a ecrire dans la base
%    resume      =  texte du resume, a ecrire dans la base
%    com         =  commentaire de l'occurence
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 17/01/2005.
% 
% 
% liste des modifications : 
%  * 17/01/2005 -> remplace rm par rm -f
%
%--------------------------------------------------------------
%
function [equi,courant,transport,generic,polar,ece,mse,resume,scronpar,com] = ...
         gcronos(numchoc,occurrence,option_str,fichier)

% initialisation des structures de cronos
data  		= [];
param 		= [];
post  		= [];
equi			= [];
courant 		= [];
transport 	= [];
generic 		= [];
polar 		= [];
ece 			= [];
mse 			= [];
resume 		= '';
scronpar    = '';
com	 		= [];		
% test des entrees
if nargin < 1
  help('gcronos');
  error('il faut donner un numero de choc !');
elseif isempty(numchoc)
  help('gcronos');
  error('il faut donner un numero de choc !');
else 
  numchoc = fix(numchoc);
end
if nargin < 2
  occurrence = 0;
end

% numero du choc avec occurence
numchoc = numchoc + rem(fix(occurrence),10)/10;

noexec = 0;
write  = 0;
if nargin < 3
   option_str ='';	
elseif strcmp(lower(option_str),'noexec')
   noexec = 1;
   option_str ='';	
elseif strcmp(lower(option_str),'write')
   write = 1;
   noexec = 1;
   option_str ='';	
end	
if nargin < 4
   fichier ='';	
end	

% donnees vides pour la sortie
% 1 - equilibre
equi.temps             = [];                 % vecteur temps associe
equi.x                 = [];                 % vecteur espace associe
equi.gcronvpr          = [];                 % element de volume [m^2]
equi.gcronspr          = [];                 % element de surface [m]
equi.gcrona            = [];                 % petit rayon equatorial [m]
equi.gcrond            = [];                 % decentrement de Shafranov [m]
equi.gcrone            = [];                 % elongation [su]
equi.scronrhom         = [];                 % rhomax [m]
equi.scrondrhomdt      = [];                 % d rhomax/ dt [m/s]
equi.gcrongrho2        = [];                 % data.equi.grho2 [su]
equi.gcrongrho         = [];                 % data.equi.grho [su]
equi.gcronr2i          = [];                 % data.equi.r2i [m^-2]
equi.gcrongrho2r2      = [];                 % data.equi.grho2r2 [m^-2]
equi.gcronri           = [];                 % data.equi.ri [m^-1]
equi.gcronrmoy         = [];                 % <R> [m]
equi.gcronfdia         = [];                 %  fonction diamagnetique [Tm]
equi.gcrongrho2b2      = [];                 % data.equi.grho2b2 [m^-2]
equi.gcronb2           = [];                 %  data.equi.b2 [T^2]
equi.gcronb2i          = [];                 %  data.equi.b2i [T^-2]
equi.gcronftrap        = [];                 %  fraction de pieges [su]
equi.scronfail         = [];                 % indicateur de convergence de l'equilibre (0 = ok)
equi.certification     = -2;                 % certification associee a la structure equilibre

% 2 - courant
courant.temps          = [];                 % vecteur temps associe
courant.x              = [];                 % vecteur espace associe
courant.gcronpsi       = [];                 % flux poloidal [Wb]
courant.gcronepar      = [];                 % champ electrique // [Vm^-1]
courant.gcronjmoy      = [];                 % courant toroidal moyen [Am^-2]
courant.gcronjeff      = [];                 % courant effectif moyen (//B) [Am^-2]
courant.gcronptot      = [];                 % pression isotropique totale [Pa]
courant.gcronq         = [];                 % facteur de securite [su]
courant.gcronjfci      = [];                 % source de courant due a FCI  [Am^-2]
courant.gcronjfce      = [];                 % source de courant due a FCE  [Am^-2]
courant.gcronjhyb      = [];                 % source de courant due a l'hybride  [Am^-2]
courant.gcronjboot     = [];                 % source de courant due au bootstrap  [Am^-2]
courant.gcroneta       = [];                 % resistivite  [Ohm m]
courant.gcronbpol      = [];                 % champ poloidal moyen   [T]
courant.scronip        = [];                 % courant  plasma   [A]
courant.certification  = -2;                 % certification associee a la structure courant

% 3 - chaleur electronique
transport.temps        = [];                 % vecteur temps associe
transport.x            = [];                 % vecteur espace associe
transport.gcronpe      = [];                 % pression electronique [Pa]
transport.gcronqe      = [];                 % flux de chaleur electronique [W m^-2]
transport.gcronse      = [];                 % source de chaleur electronique [W m^-3]
transport.gcronsefci   = [];                 % source de chaleur electronique due a FCI [W m^-3]
transport.gcronsefce   = [];                 % source de chaleur electronique due a FCE [W m^-3]
transport.gcronsehyb   = [];                 % source de chaleur electronique due a l'hybride [W m^-3]
transport.gcronseohm   = [];                 % source de chaleur ohmique [W m^-3]
transport.gcronserad   = [];                 % pertes par rayonnement [W m^-3]
transport.gcronsei     = [];                 % terme d'equipartition [W m^-3]
transport.certification= -2;                 % certification associee a la structure courant transport

% 4 - generic
generic.temps          = [];                 % vecteur temps associe
generic.x              = [];                 % vecteur espace associe
generic.gcronne        = [];                 % densite electronique [m^-3]
generic.gcronae        = [];                 % densite ionique sur densite electronique [su]
generic.gcronpion      = [];                 % pression ionique [Pa]
generic.gcronzeff      = [];                 % profil de zeff [Pa]
generic.gcronqeneo     = [];                 % flux de chaleur electronique neoclassique [W m^-2]
generic.gcronqineo     = [];                 % flux de chaleur ionique neoclassique [W m^-2]
generic.gcronsifci     = [];                 % source de chaleur ionique due a FCI [W m^-3]

% attention ce groupe a pas la meme nombre de voies
% 1  = FCE
% 2  = FCI
% 3  = HYB
% 4  = IMPUR
generic.gcronmode      = [];                 % indicateur des temps effectif pour le calcul des source et du zeff 

generic.scronnbconv    = [];                 % nombre de boucles de convergences par pas de temps pour la resolution des equations de diffusion
generic.certification  = -2;                 % certification associee a la structure de donnees generiques

% le resume
resume                 = 'Pas de donnees produites';
% struture param et mode sous forme de texte matlab
scronpar       		  = '';                 % structure param sous forme de texte

% 5 - polarimetrie
polar.temps          = [];                 % vecteur temps associe
polar.x              = [];                 % deuxieme coordonnee du groupe
polar.gcronalpha     = [];                 % angle recalculer de la polarimetrie [radian]
polar.certification  = -2;                 % certification associee 

% 6 - mse
mse.temps          = [];                 % vecteur temps associe
mse.x              = [];                 % deuxieme coordonnee du groupe
mse.gcronmse       = [];                 % angle recalculer de la mse [radian]
mse.certification  = -2;                 % certification associee 

% 7 - ECE
ece.temps          = [];                 % vecteur temps associe
ece.x              = [];                 % deuxieme coordonnee du groupe
ece.gcronrece      = [];                 % rayon de mesure de l'ece recalculer avec l'equilibre cronos [m]
ece.certification  = -2;                 % certification associee 

%
% selon le mode
%
if isempty(fichier)  & isempty(option_str)
   % execution automaitique apres tprof
	fprintf('Prepartion d''une simulation Cronos Choc T-S  #%g  (diffusion du courant seule)\n',numchoc);

	% 1 - lecture de la base temps tprof
	[ip,times] = tsbase(numchoc,'sprofip');
	if isempty(times)
		disp('Pas de donnees tprof pour ce choc ...');
		return
	end

   % lecture du scenario
	[wdia,times1] = tsbase(numchoc,'sprofwdia');
	[nl,times2] = tsbase(numchoc,'sprofnl');
	[plh,times3] = tsbase(numchoc,'sprofplh');
	 plh(plh < 0.1 ) = 0;
	 seuilhyb  = 0.3;
	 if any(plh > seuilhyb)
	    hybride = 1;
	 else 
	    hybride = 0;
	 end
	 
	[pfci,times4] = tsbase(numchoc,'gprofpfci');
	 pfci(pfci < 0.1 ) = 0;
	 if any(sum(pfci,2) > 0.3)
	    fci = 1;
	 else 
	    fci = 0;
	 end
	 
	% pour tous les chocs
         occ_ecrh   = tsoccur('tfce',fix(numchoc));
         if ~isempty(occ_ecrh )
            [pecrh,tecrh] = tsbase(fix(numchoc) + min(occ_ecrh) / 10,'gpfce');
            pecrh = interp1(tecrh,pecrh,times);
            tecrh = times;
         else
	    [pecrh,polarfce,phi_pol,phi_tor] = zecrh(fix(numchoc),times);
         end
	 if any(sum(pecrh,2) > 0.1)
	    ecrh = 1;
	 else 
	    ecrh = 0;
	 end
	
	% te0 pour la base temps
	[te,tte] = tsbase(numchoc,'gproftefit');
         te0      = mean(te(:,1:3),2);
	
	% selection de la base temps pertinente pour cronos
	if ~any(ip > 0.3)
	   disp('courant trop petit')
		return
	end
	% point de depart de la simulation
	dipdt	= pdederive(times,ip,2,2,1,1);
	if times(1) > 0.5 
	   kmin = 1;
	else
	   ind1  = find( (ip > 0.3) & (dipdt < 1) & (times < 0.75));
		if isempty(ind1)
	      kmin  = min(find( (ip > 0.3) & (dipdt < 1)));
			if isempty(kmin)
				kmin = 3;
			end
		else
			kmin  = max(ind1(find(ip(ind1) == min(ip(ind1)))));
		end
   end
	% point final de la simulation
	kmax    = max(find( (times > times(kmin)) & (ip > 0.3) & (dipdt > -2)));
	if isempty(kmax)
	   kmax = length(times);
	end
	
	% creation de la sous base temps tenant compte des evolution du choc
	% base temps initiale (120 temps max)
	dk = max(ceil((kmax-kmin)/ 120),1);
	t00    =  times(kmin:dk:kmax);
	lim    =  length(t00);
	tip    =  zextraitnoeudlim(times(kmin:kmax),ip(kmin:kmax),31,2,lim);
	idt    =  find((times2 >= min(tip)) & (times2 <= max(tip)));
	tnl    =  zextraitnoeudlim(times2(idt),nl(idt),31,2,lim);
	idt    =  find((tte >= min(tip)) & (tte <= max(tip)));
	tte0   =  zextraitnoeudlim(tte(idt),te0(idt),31,2,lim);
	idt    =  find((times1 >= min(tip)) & (times1 <= max(tip)));
	twdia  =  zextraitnoeudlim(times1(idt),wdia(idt),31,2,lim);
	if hybride == 1 
	    	idt    =  find((times3 >= min(tip)) & (times3 <= max(tip)));
		tplh     =  zextraitnoeudlim(times3(idt),plh(idt) .* (plh(idt) > seuilhyb),31,2,lim);
		tplh([1,end]) = [];
	else
	    tplh =[];	 
	end
	if fci == 1	 
	    idt    =  find((times4 >= min(tip)) & (times4 <= max(tip)));
	    tpfci  =  zextraitnoeudlim(times4(idt),sum(pfci(idt,:),2),31,2,lim);
	    tpfci([1,end]) = [];
	else
	    tpfci  = []; 
	end
	if ecrh == 1
	    tpecrh =  zextraitnoeudlim(times(kmin:kmax),sum(pecrh(kmin:kmax,:),2),31,2,lim);
		 tpecrh([1,end]) = [];
         else
	    tpecrh =  [];
	end
	
	% vecteur temps pour cronos
	temps = union(tip,t00);
	temps = union(temps,tnl);
	temps = union(temps,tte0);
	temps = union(temps,twdia);
	temps = union(temps,tplh);
	temps = union(temps,tpfci);
	temps = union(temps,tpecrh);

	% securite sur le vecteur de temps
	indbad = find(diff(temps) <=1e-4);
	if ~isempty(indbad)
		disp('negtive time step will be corrected');
		while(~isempty(indbad))
			temps(indbad) =[];
			indbad = find(diff(temps) <=1e-4);
		end	
	end
	
	
	% acces au donnees pour cronos
	info    = ztsacces;
	option  = info.valeur;
	
	% choix des option
	% le profil de courant de depart 
	option.fjli = 1;   % donnees Ip et li
	% choix du profil de ti
	option.tiop = 0;   % mesures si diponibles, sinon loi d'echelle
	% calcul neoclassique optimiser
	option.fast = 2;
	option.plotonoff   = 0;       % pas de plot
	option.verbose     = 0;       % silencieux
	option.rebuilt     = 0;       % pas reconstruction 
	option.post        = 1;       % calcul du postprocessing

	% pour 30067
	%option.faczeff = 1.5;
	
	% creation des donnees
      user_real = getenv('USER');
   if strcmp(user_real,'cgc')
	   [cr,data,param]=ztsacces(numchoc,'/usr/drfc/zineb/data',temps,option);
   else
      lieu_reel  = fullfile(getenv('HOME'),'zineb/data');
	   [cr,data,param]=ztsacces(numchoc,lieu_reel,temps,option);
   end

   % supression des sauvegardes 
	if write == 0
   	param.gene.nbsauve = inf;
   	param.gene.file    = '';
		param.gene.origine ='';
		param.gene.rapsauve = '';

   	% mode muet
		param.gene.verbose = 0;
	end
	
	% regalge pour le zeff
	param.cons.impur.zeff     =  1;
	param.cons.impur.exposant =  0;
	data.mode.zeff            =  2 .* ones(size(data.mode.zeff));
	
	% choix des modules pour les calculs de sources
	% fci
	param.fonction.fci  = 'zfcifile';
	param.cons.fci      =  zrepicparam(param.cons.fci,param.fonction.fci,param.nombre.fci);
	data.mode.fci       =  0 .* data.mode.fci;
	if fci == 1
	       %data.mode.fci(temps > min(tpfci)) = 3;
	       %data.mode.fci(temps == min(tpfci)) = 2;
			 %indfci  = interp1(temps,1:length(temps),tpfci,'nearest');
			 %data.mode.fci(indfci) = 2;
			 pfci  = interp1(times,sum(pfci,2),temps);
			 data.mode.fci(pfci > 0.3) = 2;
   end
	
	% hybride
	if isempty(strmatch('xdur',param.from.source.desc,'exact'))
		param.fonction.hyb     = 'zdelphe';
		param.cons.hyb         =  zrepicparam(param.cons.hyb,param.fonction.hyb,param.nombre.hyb);
		param.cons.hyb.conf    = 'MIXTE';
		param.cons.hyb.effmult = 0.6;
		param.cons.hyb.D       = 0;
		data.mode.hyb       =  0 .* ones(size(data.mode.hyb));
		if hybride == 1
				% le nombre d'appel a delphine doit etre limiter a 9
				plhs      = medfilt1(sum(abs(data.cons.hyb),2),5);
				lhnok     = (all(abs(data.cons.hyb)< 2e5) | all(angle(data.cons.hyb)<=1));
				plhs(lhnok) = 0; 
	    			tplh2     =  zextraitnoeudlim(data.gene.temps,plhs,31,2,11);
		 		tplh2([1,end]) = [];
	     	   		data.mode.hyb(temps > min(tplh2)) = 3;
	         		data.mode.hyb(temps == min(tplh2)) = 2;
			 	indhyb  = interp1(temps,1:length(temps),tplh2,'nearest');
			 	data.mode.hyb(indhyb) = 2;
 	  end
	else
		data.mode.hyb       =  2 .* ones(size(data.mode.hyb));
		param.cons.hyb.xdur =  1;
	end
	
	% ECRH
	param.fonction.fce  = 'zremafile';
	param.cons.fce      =  zrepicparam(param.cons.fce,param.fonction.fce,param.nombre.fce);
	data.mode.fce       =  0 .* data.mode.fce;
	if ecrh == 1
			 pfce  = interp1(times,sum(pecrh,2),temps);
			 data.mode.fce(pfce > 0.1) = 2;
   end
	
	% pas de temps d'execution
	param.split.dtmax  = 5e-2;
	param.split.dtmin  = 0.1e-3;
	param.split.equi   = 5e-2;
	
	
	% securite pour le temps initial
	kmin = min(find(all( (data.prof.pe > 0) & (data.prof.pion > 0) & (data.prof.ne > 0),2)));
	if kmin > 1
		fprintf('Probleme de profils initiaux, temps initial deplace �%g s\n',data.gene.temps(kmin));
		param.gene.kmin = kmin;
		param.gene.k    = kmin;
	end
	kmin = min(find( (data.prof.te(:,1) > 199) & (data.prof.ti(:,1) > 100) & (data.prof.ne(:,1) > 3e18)));
	if kmin > 1
		fprintf('Probleme de profils initiaux, temps initial deplace �%g s\n',data.gene.temps(kmin));
		param.gene.kmin = kmin;
		param.gene.k    = kmin;
	end
	
	
	% securite sur jmoy ini
	jforme    = (1 - param.gene.x .^ 2) .^ 2;
	iforme    = trapz( param.gene.x,jforme .* (2*pi));
	kforme    = param.gene.k;
	j0ini     = data.cons.ip(kforme) ./ data.geo.a(kforme) ./ iforme;
	jmoyini   = j0ini .* jforme;
	data.prof.jmoy(param.gene.k,:) = jmoyini;
	data.prof.jmoy(param.gene.kmin,:) = jmoyini;
												
	
	% commentaire : on uitilise les informations de bile
	% param.from.creation.com = 'Traitement tcronos';
	try 
			param.from.creation.com = sprintf('tprof: %s',param.from.shot.info.bile.comment.l1);
	catch
	      param.from.creation.com = '';
	end

	% execution de cronos
	if noexec == 0
		[cr,data_ok,data,param,post] = exec_cronos(data,param,post);
		if cr ~=  0 
			disp('Erreur d''execution de Cronos')
			resume_mem = resume;
			resume  = lasterr;
			if isempty(resume)
			   resume  = resume_mem;
			end
	 	end 
		if data_ok ==  0 
			disp('Pas de donnees produite')
			resume_mem = resume;
			resume  = lasterr;
			if isempty(resume)
			   resume  = resume_mem;
			end
			return
	 	end 
	elseif write == 1
		% ecriture des donnnes le workspace
      equi = param;
      courant = data;
      transport = post;
      generic = 0;   % cr
      %
		%zassignin('base','data',data);
		%zassignin('base','param',param);
		%zassignin('base','post',post);
		%evalin('base','cr = 0;');
		%zuirename
		return
	   %param.gene.k = param.gene.kmax;
	else 
	   % juste pour les tests
	   param.gene.k = param.gene.kmax;
	end		 

elseif isempty(fichier) 
   
	switch option_str
	case   'batch'
	       % execution batch a partir d'un fichier preparer sur le user cgc
			 % test de l'acces a cgc_data
			 [s,t] = unix('ls /usr/drfc/cgc/cgc_data/zineb/data');
			 if s ~= 0
			     % il faut monter le disque
				  ! mount -t nfs deneb:/usr/deneb/cgc/cgc_data_mnt /usr/drfc/cgc/cgc_data
				  ! sync
			 end
			 % recherche de la listes des fichiers a executer
			 [s,t] = unix(sprintf('ls -t1 /usr/drfc/cgc/cgc_data/zineb/data/TSAFAIRE/*#%d#.mat*',numchoc));
			 if s ~= 0 
			      disp('Pas de fichier a traiter');
					return
			 end
			 % il faut deplacer le fichier
			 tt    = tseparec(t);
			 file1 = deblank(tt(1,:));
			 
			 % chargement des donnees
			 [data,param,post]=zuiload(file1);
			 % suppression du fichier cree par la fonction de batch
			 [s,t] = unix(sprintf(' rm -f %s',file1));
			 if s ~= 0
			    disp('Probleme de suppression du fichier cree par la fonction de batch :');
			    disp(t)
			 end
			 % mode muet
			 param.gene.verbose = 0;
			 param.plot.onoff = 0;
			 param.gene.nbsauve = inf;

			 % execution de cronos 
			 [cr,data,param,post] = exec_cronos(data,param,post);
			 if cr ~=  0 
			     disp('Erreur d''execution de Cronos')
			 end 
			 
			 % sauvegarde du resultat
			 disp('Sauvegarde finale\n');
			% sauvegarde finale
			% compactage des donnees
			data=zreduit(param,data,'compact');
			% sauvegarde
			if  verLessThan('matlab','7.0')
						savedb(param.gene.file,'param','data','post','-V6');
			else
						savedb(param.gene.file,'param','data','post');
			end
			%save(param.gene.file,'param','data','post');
			% compression du fichier
			zgzip(param.gene.file,'compress');
			% decompactage des donnees
			data=zreduit(param,data,'uncompact');
			if ~isempty(param.gene.rapsauve)
						cr =zrm_rapsauve(param.gene.rapsauve);
						disp('Probleme lors de la suppression des fichiers intermediaires\n');
			end
	otherwise
	    disp('Option non implementee !')
		 return
	end
   
   % ce mode n'ecrit pas dans la base
	return
		
	
else
   
	% chargement du fichier
	if ischar(fichier)
        [data,param,post]=zuiload(fichier);
   elseif ~isempty(getenv('DISPLAY'))
	     % cas interractif pour les tests
        [data,param,post]=zuiload;
	end
   if isempty(data)
	    disp('Pas de donnees associees a ce fichier');
		 return
   end
	
	%
	% test du numero du choc
	%
	if  (fix(param.from.shot.num) ~= fix(numchoc) ) & (numchoc >= 200)
	   disp('Le numero du choc associe aux donnees contenues dans le fichier ne correspond pas au numero du choc demande !');
		return
	end
	
end

% extraction des indices valides
kmin = param.gene.kmin;
kmax = min(param.gene.k, param.gene.kmax);


% remplissage desstructures de donnees
% donnees vides pour la sortie
% 1 - equilibre
equi.temps             =  znonan(data.gene.temps,kmin,kmax,'temps');                 % vecteur temps associe
equi.x                 =  znonan(param.gene.x,'x');                                  % vecteur espace associe
equi.gcronvpr          =  znonan(data.equi.vpr,kmin,kmax,'vpr');                   % element de volume [m^2]
equi.gcronspr          =  znonan(data.equi.spr,kmin,kmax,'spr');                   % element de surface [m]
equi.gcrona            =  znonan(data.equi.a,kmin,kmax,'a');                     % petit rayon equatorial [m]
equi.gcrond            =  znonan(data.equi.d,kmin,kmax,'d');                     % decentrement de Shafranov [m]
equi.gcrone            =  znonan(data.equi.e,kmin,kmax,'e');                     % elongation [su]
equi.scronrhom         =  znonan(data.equi.rhomax,kmin,kmax,'rhomax');                % rhomax [m]
equi.scrondrhomdt      =  znonan(data.equi.drhomaxdt,kmin,kmax,'drhomaxdt');             % d rhomax/ dt [m/s]
equi.gcrongrho2        =  znonan(data.equi.grho2,kmin,kmax,'grho2');                 % data.equi.grho2 [su]
equi.gcrongrho         =  znonan(data.equi.grho,kmin,kmax,'grho');                  % data.equi.grho [su]
equi.gcronr2i          =  znonan(data.equi.r2i,kmin,kmax,'r2i');                   % data.equi.r2i [m^-2]
equi.gcrongrho2r2      =  znonan(data.equi.grho2r2,kmin,kmax,'grho2r2');               % data.equi.grho2r2 [m^-2]
equi.gcronri           =  znonan(data.equi.ri,kmin,kmax,'ri');                    % data.equi.ri [m^-1]
equi.gcronrmoy         =  znonan(data.equi.rmoy,kmin,kmax,'rmoy');                  % <R> [m]
equi.gcronfdia         =  znonan(data.equi.F,kmin,kmax,'F');                     %  fonction diamagnetique [Tm]
equi.gcrongrho2b2      =  znonan(data.equi.grho2b2,kmin,kmax,'grho2b2');                     % data.equi.grho2b2 [m^-2]
equi.gcronb2           =  znonan(data.equi.b2,kmin,kmax,'B^2');                     %  data.equi.b2 [T^2]
equi.gcronb2i          =  znonan(data.equi.b2i,kmin,kmax,'1/B^2');                     %  data.equi.b2i [T^-2]
equi.gcronftrap        =  znonan(data.equi.ftrap,kmin,kmax,'ftrap');                     %  fraction de pieges [su]
equi.scronfail         =  znonan(data.equi.fail + data.equi.oscil * 10,kmin,kmax,'fail');                 % indicateur de convergence de l'equilibre (0 = ok)

if isempty(equi.scronrhom)
	 equi.temps              =  [];                % vecteur temps associe
    equi.x                  =  [];                % vecteur espace associe
    equi.certification      = -2;                 % certification associee a la structure equilibre
elseif any(equi.scronfail)
    equi.certification      = -1;                 % certification associee a la structure equilibre
else
    equi.certification      = 0;                 % certification associee a la structure equilibre
end

% 2 - courant
courant.temps          = znonan(data.gene.temps,kmin,kmax,'temps');                 % vecteur temps associe
courant.x              = znonan(param.gene.x,'x');                              % vecteur espace associe
courant.gcronpsi       = znonan(data.prof.psi,kmin,kmax,'psi');                   % flux poloidal [Wb]
courant.gcronepar      = znonan(data.prof.epar,kmin,kmax,'epar');                  % champ electrique // [Vm^-1]
courant.gcronjmoy      = znonan(data.prof.jmoy,kmin,kmax,'jmoy');                  % courant toroidal moyen [Am^-2]
courant.gcronjeff      = znonan(data.prof.jeff,kmin,kmax,'jeff');                  % courant effectif moyen (//B) [Am^-2]
courant.gcronptot      = znonan(data.prof.ptot,kmin,kmax,'ptot');                  % pression isotropique totale [Pa]
courant.gcronq         = znonan(data.prof.q,kmin,kmax,'q');                     % facteur de securite [su]
courant.gcronjfci      = znonan(data.source.fci.j,kmin,kmax,'jfci');               % source de courant due a FCI  [Am^-2]
courant.gcronjfce      = znonan(data.source.fce.j,kmin,kmax,'jfce');               % source de courant due a FCE  [Am^-2]
courant.gcronjhyb      = znonan(data.source.hyb.j,kmin,kmax,'jhyb');               % source de courant due a l'hybride  [Am^-2]
courant.gcronjboot     = znonan(data.source.jboot,kmin,kmax,'jboot');               % source de courant due au bootstrap  [Am^-2]
courant.gcroneta       = znonan(data.coef.eta,kmin,kmax,'eta');                   % resistivite  [Ohm m]
courant.gcronbpol      = znonan(data.prof.bpol,kmin,kmax,'bpol');                  % champ poloidal moyen   [T]
courant.scronip        = znonan(data.gene.ip,kmin,kmax,'ip');                    % courant  plasma   [A]

if isempty(courant.gcronpsi) 
	courant.temps          = [];                 % vecteur temps associe
	courant.x              = [];                              % vecteur espace associe
	courant.certification  = -2;                 % certification associee a la structure courant
elseif (sum(data.gene.conv)/ param.gene.nbt) > 1e2
	courant.certification  = -1;                 % certification associee a la structure courant
else
	courant.certification  = 0;                 % certification associee a la structure courant
end

% 3 - chaleur electronique
transport.temps        = znonan(data.gene.temps,kmin,kmax,'temps');                 % vecteur temps associe
transport.x            = znonan(param.gene.x,'x');                              % vecteur espace associe
transport.gcronpe      = znonan(data.prof.pe,kmin,kmax,'pe');                    % pression electronique [Pa]
transport.gcronqe      = znonan(data.prof.flux.qe,kmin,kmax,'qe');               % flux de chaleur electronique [W m^-2]
transport.gcronse      = znonan(data.source.totale.el,kmin,kmax,'sel');           % source de chaleur electronique [W m^-3]
transport.gcronsefci   = znonan(data.source.fci.el,kmin,kmax,'sefci');              % source de chaleur electronique due a FCI [W m^-3]
transport.gcronsefce   = znonan(data.source.fce.el,kmin,kmax,'sefce');              % source de chaleur electronique due a FCE [W m^-3]
transport.gcronsehyb   = znonan(data.source.hyb.el,kmin,kmax,'sehyb');              % source de chaleur electronique due a l'hybride [W m^-3]
transport.gcronseohm   = znonan(data.source.ohm,kmin,kmax,'ohm');                 % source de chaleur ohmique [W m^-3]
transport.gcronserad   = znonan(data.source.prad + data.source.brem + data.source.cyclo,kmin,kmax,'prad');                 % pertes par rayonnement [W m^-3]
transport.gcronsei     = znonan(data.source.qei  + data.source.qneo,kmin,kmax,'qei');                 % terme d'equipartition [W m^-3]

if isempty(transport.gcronpe) 
	transport.temps         = [];                 % vecteur temps associe
	transport.x             = [];                 % vecteur espace associe
	transport.certification = -2;                 % certification associee a la structure courant transport
elseif (param.from.option.tiop == 0) | (param.from.option.tiop == 5)
   % tiop = 0 -> tiprof
   % tiop = 5 -> fit des mesure par tprof
	 transport.certification = 0;                 % certification associee a la structure courant transport
else
	 transport.certification = -1;                 % certification associee a la structure courant transport
end
% 4 - generic
generic.temps          = znonan(data.gene.temps,kmin,kmax,'temps');                 % vecteur temps associe
generic.x              = znonan(param.gene.x,'x');                     % vecteur espace associe
generic.gcronne        = znonan(data.prof.ne,kmin,kmax,'ne');                    % densite electronique [m^-3]
generic.gcronae        = znonan(data.prof.ae,kmin,kmax,'ae');                    % densite ionique sur densite electronique [su]
generic.gcronpion      = znonan(data.prof.pion,kmin,kmax,'pion');                  % pression ionique [Pa]
generic.gcronzeff      = znonan(data.prof.zeff,kmin,kmax,'zeff');                  % profil de zeff [Pa]
generic.gcronqeneo     = znonan(data.neo.flux.qe,kmin,kmax,'neoqe');                % flux de chaleur electronique neoclassique [W m^-2]
generic.gcronqineo     = znonan(data.neo.flux.qion,kmin,kmax,'neoqi');              % flux de chaleur ionique neoclassique [W m^-2]
generic.gcronsifci     = znonan(data.source.fci.ion,kmin,kmax,'sifci');             % source de chaleur ionique due a FCI [W m^-3]

% attention ce groupe a pas la meme nombre de voies
% 1  = FCE
% 2  = FCI
% 3  = HYB
% 4  = IMPUR
generic.gcronmode      = znonan(cat(2,data.mode.fce,data.mode.fci,data.mode.hyb,data.mode.impur),kmin,kmax,'mode');  % indicateur des temps effectif pour le calcul des source et du zeff 

generic.scronnbconv    = znonan(data.gene.conv,kmin,kmax,'conv');        % nombre de boucles de convergences par pas de temps pour la resolution des equations de diffusion

divers.source.fci.err  = data.source.fci.err;
divers.post.polar = [];
if isfield(post,'polar')
   if isfield(post.polar,'nl')
      for klk = 1:size(post.polar.nl,2)
         divers.post.polar = setfield(divers.post.polar,sprintf('nl%d',klk),znonan(post.polar.nl(:,klk),kmin,kmax,'nl'));
      end
   end
end

scronpar               = zpar2str(param,data.mode,data.cons,divers,kmin,kmax);     % structure param sous forme de texte

if isempty(generic.temps)
	generic.temps          = [];                 % vecteur temps associe
	generic.x              = [];                 % vecteur espace associe
	generic.certification  = -2;                 % certification associee a la structure de donnees generiques
else
	generic.certification  = 0;                 % certification associee a la structure de donnees generiques
end

% commentaire de l'occurrence
com = sprintf('%s:	%s',param.from.creation.user,param.from.creation.com);
if isempty(param.from.creation.com)
   try 
			com = sprintf('tprof: %s',param.from.shot.info.bile.comment.l1);
	catch
	      com ='';
	end
else
	com = sprintf('%s:	%s',param.from.creation.user,param.from.creation.com);
end
if isempty(com)
	com ='execution automatique de Tcronos'; 
end

% 5 - polarimetrie
if isfield(post,'polar')
     if size(post.polar.af,1) == length(data.gene.temps)
           polar.temps          = znonan(data.gene.temps,kmin,kmax,'temps');               % vecteur temps associe
           polar.gcronalpha     = znonan(post.polar.af,kmin,kmax,'af');                 % angle recalculer de la polarimetrie [radian]
           polar.x              = 1:size(polar.gcronalpha,2);               % indices
           polar.certification  = courant.certification;                 % certification associee 
	  else
          polar.certification  = -2;                 % certification associee 
	  end
else
     polar.certification  = -2;                 % certification associee 
end

% 6 - mse
% resever pour plus tard
mse.certification  = -2;                 % certification associee 

% 7 - ECE
if isfield(post,'ece')
     if size(post.ece.Rece,1) == length(data.gene.temps)
           ece.temps            = znonan(data.gene.temps,kmin,kmax,'temps');               % vecteur temps associe
           ece.gcronrece        = znonan(post.ece.Rece,kmin,kmax,'rece');                  % rayon de mesure de l'ece recalculer avec l'equilibre cronos [m]
           ece.x                = 1:size(ece.gcronrece,2);   % indices
           ece.certification  = equi.certification;                    % certification associee 
	  else
          ece.certification  = -2;                 % certification associee 
	  end
else
     ece.certification  = -2;                 % certification associee 
end

% 8 - creation du resume
[cr,message] = zresume(param,data,post,'zou');
if cr ~= 0
   resume = '';
else
   resume = message;
   % ajout des commentaire de tprof
	try
		resume = sprintf('%s\n---------------------------------------------------------\n',resume);
		resume = sprintf('%sCommentaires en provenance de tprof :\n',resume);
		if isfield(param.from.shot.info.bile.comment,'l1')
			resume = sprintf('%s%s\n',resume,param.from.shot.info.bile.comment.l1);
		end 
		if isfield(param.from.shot.info.bile.comment,'l2')
			resume = sprintf('%s%s\n',resume,param.from.shot.info.bile.comment.l2);
		end 
		if isfield(param.from.shot.info.bile.comment,'l3')
			resume = sprintf('%s%s\n',resume,param.from.shot.info.bile.comment.l3);
		end 
		if isfield(param.from.shot.info.bile.comment,'l4')
			resume = sprintf('%s%s\n',resume,param.from.shot.info.bile.comment.l4);
		end 
		if isfield(param.from.shot.info.bile.comment,'l5')
			resume = sprintf('%s%s\n',resume,param.from.shot.info.bile.comment.l5);
		end 
	end
end


% fin du traitement

% fonction de protecttion contre les nans et Inf
function s = znonan(e,k1,k2,name)

% extraction de la sous matrice
if nargin >= 3
   e = e(k1:k2,:);
elseif nargin == 2
  name = k1;
end

ind = find(~ isfinite(e));
if ~isempty(ind)
	e(ind) = zeros(1,length(ind));
   fprintf('Attention : supression de %d NaN ou Inf dans %s\n',length(ind),name);
end

ind = find(imag(e));
if ~isempty(ind)
	e(ind) = real(e);
   fprintf('Attention : supression de %d Imag dans %s\n',length(ind),name);
end

s = e;


function  [xx,yy]= zextraitnoeudlim(x,y,nb,mode,lim)

if nargin < 5
 lim = inf;
end

% anti retour an arriere
if (mode == 0) | (mode == 2)
	indnok = find(diff(x)<=0);
	while(~isempty(indnok))
		x(indnok) = [];
		y(indnok) = [];
		indnok = find(diff(x)<=0);
	end
end

[xx,yy]= zextraitnoeud(x,y,nb,mode);
if mode == 1 
	return
elseif ~isfinite(lim)
	return
else
   while (length(xx) > lim) & (nb > 1)
	    if nb == 2
		 	nb = 1;
		 else 
	    	nb = fix(nb / 2) + 1;
		 end
		 [xx,yy]= zextraitnoeud(x,y,nb,mode);
   end

end
