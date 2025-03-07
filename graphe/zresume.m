%
% cree la page de resume de cronos
% cr = zresume(param,data,post);
% [cr,message] = zresume(param,data,post,zou);
%
function [cr,message] =zresume(p1,p2,p3,zou)

% initialisation de cr
cr =0;

% gestion des entrees
file =[];
edit =0;
efface = 0;
if nargin >= 2
   % mode affichage seul
   param = p1;
   data  = p2; 
   if nargin == 3
      post =p3;
   else
      post =[];
   end
   % fichier temporaire
   file = tempname;
   % affichage en final
   edit =1;
   % efface le fichier apres execution
   efface =1;
elseif nargin == 0
   % mode selection du fichier + cretion fichier resume
   [data,param,post]=zuiload;
   % nom du fichier 
   file = param.gene.file;
   % affichage en final
   edit =1;
elseif isempty(p1)
   % mode selection du fichier + creation fichier resume
   [data,param,post]=zuiload;
   % nom du fichier 
   file = param.gene.file;
   % affichage en final
   edit =1;
else
   % chargement des donnees + creation du fichier resume
   [data,param,post]=zuiload(p1);
   % nom du fichier 
   file = p1;
   % variable fichier
   fichier = p1;
end

% test du contenu
if ~isstruct(param) | ~isstruct(data)
   disp('Format de fichier incorrect !')
   disp('Reprise d''execution impossible')
   cr = -10094;
   return
end

% formatage du nom du fichier
% supression du .gz
[pf,nf,ext,ver]=fileparts(file);
file = fullfile(pf,nf);
% supression du .mat
[pf,nf,ext,ver]=fileparts(file);
file = fullfile(pf,nf);
fichier = file;
% ajout extention
file = strcat(file,'.resume_cronos');

% ouverture du fichier resume
[fid,mess] = fopen(file,'w');
if fid <3
  fprintf('Erreur lors de la creation du fichier %s :\n%s\n',file,mess);
  cr = -101;
  return
end

% calcul du temps cpu
cpu =real(data.gene.cputime);
ind = find(isfinite(cpu));
if ~isempty(ind)
   if length(ind) > 1
     dcpudt = zdxdt(cpu(ind),data.gene.temps(ind));
     cpu = cumtrapz(data.gene.temps(ind),dcpudt .* (dcpudt >=0));
   else
     cpu = cpu(ind);
   end
end
if length(cpu) < length(data.gene.temps)
   cpu((end+1):length(data.gene.temps)) = cpu(end);
end
dd =real(data.gene.datation);
ind = find(isfinite(dd));
if ~isempty(ind)
   if length(ind) > 1
   ddddt = zdxdt(dd(ind),data.gene.temps(ind));
   dd = cumtrapz(data.gene.temps(ind),ddddt .* (ddddt >=0));
   else
     dd = dd(ind);
   end
end
if length(dd) < length(data.gene.temps)
   dd((end+1):length(data.gene.temps)) = dd(end);
end


% 1- titre et identite
user = strrep(param.from.creation.user(1,:),sprintf('\n'),'');
langue = getappdata(0,'langue_cronos');
if strcmp(langue,'anglais')
  fprintf(fid,'Shot %s #%.1f de %g a %g s (%s) | %s \n',param.from.machine,param.from.shot.num, ...
        param.gene.tdeb,param.gene.tfin,zdate(param.from.shot.date),param.gene.filetype);
  fprintf(fid,'Summary of a CRONOS run of %s (Ver %3.1f); %s writes :\n%s\n', ...
        zdate(param.gene.date_exec),param.gene.version_zineb,user,param.from.creation.com);
	
% 2 - base temps
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  fprintf(fid,'time slice (s) : #%d times, tmin= %.3g, tmax= %.3g, dtmin= %.3g, dtmax= %.3g\n', ...
         length(data.gene.temps),min(data.gene.temps),max(data.gene.temps), ...
	 min(diff(data.gene.temps)),max(diff(data.gene.temps)));
  fprintf(fid,'last time calculation : %.3g  (cpu = %g s, machine = %g s)\n', ...
       data.gene.temps(param.gene.k),cpu(param.gene.k) - cpu(param.gene.kmin), ...
       dd(param.gene.k) - dd(param.gene.kmin) );

% 3 - scenario
  mode = mean(data.geo.mode);
  switch mode
   case 0
    smode = 'symmetrical plasma';
   case 1 
    smode = 'up and down triangularity';
   case 2
    smode = 'separatrix';
   otherwise
    smode = 'complex';
  end
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  fprintf(fid,['Scheme : <R0>= %.3g m, <a>= %.3g m, <elong>= %.3g, <tri>= %.3g, [%s]\n', ...
             '           <Btor>= %.3g T, <Ip>= %.3g A, sign(Ip) = %d, sign(B0) = %d\n'], ...
        mean(data.geo.r0),mean(data.geo.a),mean(data.geo.e1),mean(0.5 .* (data.geo.trh1 + data.geo.trb1)), ...
	smode,mean(data.geo.b0),mean(data.cons.ip),param.gene.signe.ip,param.gene.signe.b0);
% FCE
  fprintf(fid,'ECRH (#%d) : Pmax = ',param.nombre.fce);fprintf(fid,'%4.1g ',max(abs(data.cons.fce))/1e6);	
  fprintf	(fid,'MW, Etot = %.3g MJ\n',trapz(data.gene.temps,sum(abs(data.cons.fce),2))/1e6);
% FCI
  fprintf(fid,'ICRH (#%d) : Pmax = ',param.nombre.fci);fprintf(fid,'%4.1g ',max(abs(data.cons.fci))/1e6);	
  fprintf	(fid,'MW, Etot = %.3g MJ\n',trapz(data.gene.temps,sum(abs(data.cons.fci),2))/1e6);
% Hyb
  fprintf(fid,'LH (#%d) : Pmax = ',param.nombre.hyb);fprintf(fid,'%4.1g ',max(abs(data.cons.hyb))/1e6);	
  fprintf	(fid,'MW, Etot = %.3g MJ\n',trapz(data.gene.temps,sum(abs(data.cons.hyb),2))/1e6);
% IDN
  if strcmp(param.from.machine,'JET')&(param.nombre.idn==16)
   fprintf(fid,'NBI (#%d) : Pmax = ',param.nombre.idn);fprintf(fid,'%4.1g ', ...
           [max(sum(real(data.cons.idn(:,1:8)),2))/1e6,max(sum(real(data.cons.idn(:,2:16)),2))/1e6]);	
   fprintf(fid,'MW, Etot = %.3g MJ\n',trapz(data.gene.temps,sum(real(data.cons.idn),2))/1e6);
  else
   fprintf(fid,'NBI (#%d) : Pmax = ',param.nombre.idn);fprintf(fid,'%4.1g ',max(real(data.cons.idn))/1e6);	
   fprintf(fid,'MW, Etot = %.3g MJ\n',trapz(data.gene.temps,sum(real(data.cons.idn),2))/1e6);
  end

% 4 - Gaz :
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  nb = fprintf(fid,'Gaz (Z,A) : ');
  lg = {'Maj','Min1','Min2','Imp1','Imp2'};
  for k =1:length(param.compo.z)
	if k <= length(lg)
		leg = lg{k};
	else
		leg = sprintf('Gaz%d',k);
	end
	if (param.compo.z(k) == 1) & (param.compo.a(k) ==1)
		nb = fprintf(fid,'%s = H, ',leg) + nb;
	elseif (param.compo.z(k) == 1) & (param.compo.a(k) ==2)
		nb = fprintf(fid,'%s = D, ',leg) + nb;
	elseif (param.compo.z(k) == 1) & (param.compo.a(k) ==3)
		nb = fprintf(fid,'%s = T, ',leg) + nb;
	elseif (param.compo.z(k) == 2) & (param.compo.a(k) ==3)
		nb = fprintf(fid,'%s = He3, ',leg) + nb;
	elseif (param.compo.z(k) == 2) & (param.compo.a(k) ==4)
		nb = fprintf(fid,'%s = He, ',leg) + nb;
	elseif (param.compo.z(k) == 6) & (param.compo.a(k) ==12)
		nb = fprintf(fid,'%s = C, ',leg) + nb;
	elseif (param.compo.z(k) == 8) & (param.compo.a(k) ==16)
		nb = fprintf(fid,'%s = O, ',leg) + nb;
	else
		nb = fprintf(fid,'%s = (%d,%d), ',leg,param.compo.z(k),param.compo.a(k)) + nb;
	end
	
	if (nb > 70) & (k <length(param.compo.z))
		fprintf('\n');
		nb =0;
	end
  end
%fprintf(fid,'Gaz :    Maj        Min1        Min2        Imp1        Imp2\n');
%fprintf(fid,'  Z :');fprintf(fid,'    %2d      ',param.compo.z);
%fprintf(fid,'\n  A :');fprintf(fid,'    %2d      ',param.compo.a);
  fprintf(fid,'\n<N0>:');
  for k =1:size(data.impur.impur,3)
   ss = sprintf('%8.3g',mean(data.impur.impur(:,1,k)));
   fprintf(fid,' %9.9s ',ss);
  end
  fprintf(fid,' m^-3');

  fprintf(fid,'\nZeff : <Zeffm>= %.3g, min(Zeffm)= %.3g, max(Zeffm)= %.3g\n', ...
        mean(data.cons.zeffm),min(data.cons.zeffm),max(data.cons.zeffm));
  fprintf(fid,'nH/nD: <nhnd> = %.3g, min(nhnd) = %.3g, max(nhnd) = %.3g\n', ...
        mean(data.cons.nhnd),min(data.cons.nhnd),max(data.cons.nhnd));
  fprintf(fid,'Nebar: <Nebar> = %.3g, min(Nebar) = %.3g, max(Nebar) = %.3g\t(1e19 m^-3)\n', ...
        mean(data.gene.nbar/1e19),min(data.gene.nbar/1e19),max(data.gene.nbar/1e19));

% 5- generation de courant :
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  fprintf(fid,'max of:   Ip       Icd       Ini       Iboot       ILH       Inbi       IICRH\n');
  fprintf(fid,'(MA) %6.3g   %6.3g    %6.3g    %6.3g      %6.3g     %6.3g     %6.3g\n', ...
             max(data.gene.ip)/1e6,max(data.gene.icd)/1e6,max(data.gene.ini)/1e6,max(data.gene.iboot)/1e6, ...
	     max(data.gene.ihyb)/1e6,max(data.gene.iidn)/1e6,NaN);

% 5- configuration calcul :
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  fprintf(fid,'Mode   : Psi = %s, Pe = %s, Pion = %s, Ne = %s \n',zmode(data.mode.psi,langue),zmode(data.mode.pe,langue), ...
        zmode(data.mode.pion,langue),zmode2(data.mode.nel,langue));
  spsi   = zlimite(data.mode.cons.psi,{'Ip','Vloop','Psi1','?'});
  spe    = zlimite(data.mode.cons.pe,{'Te1','qe1','Pe1','?'});
  spion  = zlimite(data.mode.cons.pion,{'Ti1','qion1','Pion1','?'});
  sne    = zlimite(data.mode.cons.ne,{'Ne1','ge1','?','?'});
  fprintf(fid,'Limit : Psi = %s, Pe = %s, Pion = %s, Ne = %s \n',spsi,spe,spion,sne);

  switch param.gene.cn
  case 0
    scn = 'implicit';
  case 0.5
    scn = 'C-N';
  case 1
    scn = 'explicit';
  otherwise
    scn = 'other';
  end
  if param.gene.adiabatic  == 1
    sadia ='adiabatic';
  else
    sadia ='explicit';
  end
  fprintf(fid,'Solver : %s, amorti = %g, init = %s, d(adia) = %g, nmax = %d\n', ...
        scn,param.gene.amorti,sadia,param.gene.delta_adia,param.gene.nmax);

  switch param.gene.force
  case 0
    sfo = 'standard';
  otherwise
    sfo = 'forced';
  end
  if param.gene.nonneg  == 0
    spos ='free';
  else
    spos ='positive';
  end
  if param.gene.guido  == 0
    gg ='free';
  else
    gg ='positive';
  end
  fprintf(fid,'Special : convergence %s, sign(Pe, Pion & Ne) -> %s, J(0) %s \n', sfo,spos,gg);

  if param.gene.modecoef == 1
    scoef = 'convectif';
  else
    scoef = 'full';
  end
  if param.gene.self == 1
    sself = 'selfconsitant';
  else
    sself = 'standard';
  end
  if param.gene.fast == 0
     sfast = 'selfconstitant';
  elseif param.gene.fast == 1
     sfast = 'fast';
  else
    sfast = 'Jdiff';
  end	 	 
  fprintf(fid,'Eq : lambda= %g, #eq= %d, coef= %s, neo= %s, sources = %s\n',param.gene.lambda, ...
        param.gene.nbeq_mode,scoef,sfast,sself);
  switch param.gene.psiequi
  case 0
   sconv ='never';
  case 1
   sconv ='always';
  case 2
   sconv ='no-convergence';
  otherwise
   sconv ='?';
  end
  fprintf(fid,'Equi : d(jmoy) = %g, amorti = %g, nmax = %d, equi -> diff = %s\n',param.gene.djmoy, ...
        param.gene.mjmoy,param.gene.nequi,sconv);

% manque sourcebord

  fprintf(fid,'Init : d(psi) = %g, nmax = %d\n',param.gene.dpsi_ini,param.gene.nequi_ini);
  fprintf(fid,'Criterium : Psi = %g, Pe = %g, Pion = %g, Ne =%g\n',param.gene.critere.psi, ...
         param.gene.critere.pe,param.gene.critere.pion,param.gene.critere.ne);


% les modules
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  nom       = fieldnames(param.fonction);
  for k = 1:size(nom,1)
	if strcmp(nom{k},'mhd')
  
	else
		% paramatre du module
		parametre = getfield(param.cons,nom{k});
		% nom du module
		fonction  = getfield(param.fonction,nom{k});
		% mode associe s'il existe
		if isfield(data.mode,nom{k})
			mode   = getfield(data.mode,nom{k});
			switch nom{k}
			case {'glacon','plot'}
				smode = zlimite(mode,{'stop','run','?','?'});
			otherwise
				smode  = zmode(mode,langue);
			end
		else
			smode   = '';
		end  
		
		% info complementaires  
		switch nom{k}
		case 'glacon'
			icomp = sprintf(', #%d',sum(data.cons.glacon(:)));
		case 'hyb'
			if isfield(parametre,'xdur')
				if getfield(parametre,'xdur') == 1
					icomp = ', xdur';
				else
					icomp = '';
				end
			end
		otherwise
			icomp = '';
		end
		
		% nomnbre de coupleur, injecteur ...
		if isfield(param.nombre,nom{k})
			nb     = getfield(param.nombre,nom{k});
		else
			nb     = 1;
		end
		if isempty(fonction)
			sparam = '';
		elseif ~isempty(parametre) 
			% parametre par defaut
			try
				info      = feval(fonction,nb);
				defaut    = info.valeur;
			catch
			   defaut    = [];    
			end
			% comparaison
			sparam = 'standard';
			pnom   = fieldnames(parametre);
			for l  = 1:size(pnom,1)
				pc  = getfield(parametre,pnom{l});
				if isfield(defaut,pnom{l})
					pd  = getfield(defaut,pnom{l});
					if length(pc) < length(pd)
					  pc(length(pd))=0;
					end
				elseif ischar(pc)
					pd = '';
				else
					pd =[];
				end
				
				if ischar(pc)
					if ~strcmp(pc,pd)
						sparam = 'customized';
						if length(icomp) < 16
							icomp = [icomp,', ',pnom{l}];
						elseif ~findstr(icomp,'...')
							icomp = [icomp,', ...'];
						end
					end
				elseif  ~isempty(pc) & ~isempty(pd) & isnumeric(pc)
					if ~all(size(pc) == size(pd))
						sparam = 'customized';
						if length(icomp) < 16
							icomp = [icomp,', ',pnom{l}];
						elseif ~findstr(icomp,'...')
							icomp = [icomp,', ...'];
						end
					elseif ~all(all(pc==pd))
						sparam = 'customized';
						if length(icomp) < 16
							icomp = [icomp,', ',pnom{l}];
						elseif ~findstr(icomp,'...')
							icomp = [icomp,', ...'];
						end
					end
				end
			end 
		else
			sparam = 'without';
		end
		
		if isempty(fonction)
			% rien
		elseif isempty(smode)
			fprintf(fid,'%8.8s : %12.12s (%12.12s)%s\n',nom{k},fonction,sparam,icomp); 
		else
			fprintf(fid,'%8.8s : %12.12s (%12.12s) -> %8.8s %s\n',nom{k},fonction,sparam,smode,icomp); 
		end
	end	 
  end

% coef de transport
  fprintf(fid,'-------------------------------------------------------------------------------\n');
% a - les multiplicateur
  fprintf(fid,'neo *    : ');
  neonom = fieldnames(param.cons.neomulti);
  for k =1:length(neonom)
	if k <length(neonom)
		fprintf(fid,'%2s = %4.2g, ',neonom{k},getfield(param.cons.neomulti,neonom{k}));
	else
		fprintf(fid,'%2s = %4.2g',neonom{k},getfield(param.cons.neomulti,neonom{k}));
	end	
  end
  fprintf(fid,'\n');
  coefnom = fieldnames(data.coef);
  smode ={};
  for k =1:length(coefnom)
	if isfield(data.mode,coefnom{k})
		mode   = getfield(data.mode,coefnom{k});
		smode{k}  = zmode(mode,langue);
		
	else
		smode{k}='';
	end
  end
  if ~isempty(strmatch('zero',smode))
	nb = fprintf(fid,'fill to zero   : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'zero')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('datas',smode))
	nb = fprintf(fid,'datas  : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'datas')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('calculation',smode))
	nb = fprintf(fid,'calculation  : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'calculation')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('complex',smode))
	nb = fprintf(fid,'complex : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'complex')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
% les modes
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  coefnom = fieldnames(data.mode);
  smode ={};
  exclusion = {'pe','pion','nel','psi','first time'};
  for k =1:length(coefnom)
	if ~isstruct(getfield(data.mode,coefnom{k}))
		if ~isfield(param.fonction,coefnom{k}) &  ...
		   ~isfield(data.coef,coefnom{k}) &  ...	
		   isempty(strmatch(coefnom{k},exclusion,'exact'))
			mode   = getfield(data.mode,coefnom{k});
			smode{k}  = zmode(mode,langue);
		else
			smode{k}='';
		end
	else
			smode{k}='';
	end
  end
  if ~isempty(strmatch('zero',smode))
	nb = fprintf(fid,'a zero   : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'zero')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('datas',smode))
	nb = fprintf(fid,'datas  : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'datas')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('calculation',smode))
	nb = fprintf(fid,'calculation  : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'calculation')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('complex',smode))
	nb = fprintf(fid,'complex : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'complex')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
% option de creation
  fprintf(fid,'-------------------------------------------------------------------------------\n');
% parametre par defaut
  if ~isfield(param.from,'creator')
    param.from.createur ='';
  end
    
  if ~isempty(param.from.createur)
	createur  = param.from.createur;
	if strcmp(createur,'zvplusacces')
	  createur = 'zplusacces';
	end
	info      = feval(createur);
	defaut    = info.valeur;
	parametre = param.from.option;
	nb = fprintf(fid,'creator : %s -> ',createur);
	% comparaison
	pnom   = fieldnames(parametre);
	for l  = 1:size(pnom,1)
		pc  = getfield(parametre,pnom{l});
		if isfield(defaut,pnom{l})
			pd  = getfield(defaut,pnom{l});
		elseif ischar(pc)
			pd = '';
		else
			pd =[];
		end
		if ischar(pc)
			if ~strcmp(pc,pd)
				nb = nb + fprintf(fid,'%s = %s, ',pnom{l},pc);
			end
		elseif isnumeric(pc)
			if  ~isempty(pc) & ~isempty(pd) & ~all(all(pc==pd))
				nb = nb + fprintf(fid,'%s = %g, ',pnom{l},pc);
			end
		end
		if nb > 70
			fprintf(fid,'\n');
			nb = 0;
		end
	end	
  else
	nb = fprintf(fid,'creator : ?');
  end 
             
             
% sortie 
%x -fichier
  fprintf(fid,'\n-------------------------------------------------------------------------------\n');
  fprintf(fid,'associated file :\t%s\n',fichier);
  [s,t]=unix(['ls ',param.gene.origine,'.*']);
  if s ~= 0
	fprintf(fid,'the file defined in param.gene.origine does not exist !\n');
  elseif strcmp(param.gene.filetype,'resultat')
	fprintf(fid,'output file : \n  %s\n',param.gene.origine);
  else
	fprintf(fid,'input file : \n  %s\n',param.gene.origine);
  end
  [s,t]=unix(['ls ',param.gene.file,'.*']);
  if s ~= 0
    fprintf(fid,'the file defined in param.gene.origine does not exist !\n');
  else
    fprintf(fid,'ouput file : \n  %s\n',param.gene.file); 
  end

% recherche du status
  if isempty(param.gene.rapsauve)
	fprintf(fid,'no fast saving\n');
  else
	rap = param.gene.rapsauve;
	ind = max(findstr(rap,'data/'));
	if ~isempty(ind)
		rap = strcat('...',rap((ind + 4):end));
	end
	
	[s,t] = unix(['ls -1 ',param.gene.rapsauve,'*']);
	if s ~= 0
		fprintf(fid,'fast saving name file (empty) : \n  %s\n',rap);
	else
		nb = size(tseparec(t),1);
		fprintf(fid,'fast saving name file (%d times to reconstruct): \n  %s\n',nb,rap);
	end
  end
% fermeture du fichier si pas sortie standart
  if fid > 2
    fclose(fid);
  end
% affichage du resume si necessaire
  if nargin  == 4
    [s,message] = unix(['cat ',file]);
    if s ~= 0
      disp('problem encountered during reading summary :')
		disp(message);
		message ='';
    end
    [s,t] = unix(['rm -f ',file]);
    if s ~= 0
        disp('Problem deleting the summary temporary file :')
		disp(t)
    end
  elseif edit == 1
     if efface == 1
       [s,t] = unix([getappdata(0,'editeur'),' ',file,' &']);
       if strcmp(computer,'ALPHA')
         [s,t] = unix(['sleep 60;rm -f ',file,'&']);
       end
     else
       [s,t] = unix([getappdata(0,'editeur'),' ',file,' &']);
     end
     if s ~= 0
       disp('unix problem :')
       disp(t)
       cr = 107;
       return
     end 
  elseif isempty(getenv('DISPLAY')) & ~isempty( param.from.creation.user)
  % envoi par mail si batch
    [s,message] = unix(['cat ',file]);
    cr_mail=zmail(strtok(param.from.creation.user),sprintf('Cronos summary : %s',file),message);
  end

else
  fprintf(fid,'Choc %s #%.1f de %g a %g s (%s) | %s \n',param.from.machine,param.from.shot.num, ...
        param.gene.tdeb,param.gene.tfin,zdate(param.from.shot.date),param.gene.filetype);
  fprintf(fid,'Resume d''une simulation Cronos du %s (Ver %3.1f); %s a dit :\n%s\n', ...
        zdate(param.gene.date_exec),param.gene.version_zineb,user,param.from.creation.com);
	
% 2 - base temps
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  fprintf(fid,'base temps (s) : #%d temps, tmin= %.3g, tmax= %.3g, dtmin= %.3g, dtmax= %.3g\n', ...
         length(data.gene.temps),min(data.gene.temps),max(data.gene.temps), ...
	 min(diff(data.gene.temps)),max(diff(data.gene.temps)));
  fprintf(fid,'temps de fin du calcul : %.3g  (cpu = %g s, machine = %g s)\n', ...
       data.gene.temps(param.gene.k),cpu(param.gene.k) - cpu(param.gene.kmin), ...
       dd(param.gene.k) - dd(param.gene.kmin) );

% 3 - scenario
  mode = mean(data.geo.mode);
  switch mode
   case 0
    smode = 'symetrique';
   case 1 
    smode = 'asymetrique';
   case 2
    smode = 'separatrice';
   otherwise
    smode = 'complexe';
  end
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  fprintf(fid,['Scenario : <R0>= %.3g m, <a>= %.3g m, <elong>= %.3g, <tri>= %.3g, [%s]\n', ...
             '           <Btor>= %.3g T, <Ip>= %.3g A, signe(Ip) = %d, signe(B0) = %d\n'], ...
        mean(data.geo.r0),mean(data.geo.a),mean(data.geo.e1),mean(0.5 .* (data.geo.trh1 + data.geo.trb1)), ...
	smode,mean(data.geo.b0),mean(data.cons.ip),param.gene.signe.ip,param.gene.signe.b0);
% FCE
  fprintf(fid,'FCE (#%d) : Pmax = ',param.nombre.fce);fprintf(fid,'%4.1g ',max(abs(data.cons.fce))/1e6);	
  fprintf	(fid,'MW, Etot = %.3g MJ\n',trapz(data.gene.temps,sum(abs(data.cons.fce),2))/1e6);
% FCI
  fprintf(fid,'FCI (#%d) : Pmax = ',param.nombre.fci);fprintf(fid,'%4.1g ',max(abs(data.cons.fci))/1e6);	
  fprintf	(fid,'MW, Etot = %.3g MJ\n',trapz(data.gene.temps,sum(abs(data.cons.fci),2))/1e6);
% Hyb
  fprintf(fid,'HYB (#%d) : Pmax = ',param.nombre.hyb);fprintf(fid,'%4.1g ',max(abs(data.cons.hyb))/1e6);	
  fprintf	(fid,'MW, Etot = %.3g MJ\n',trapz(data.gene.temps,sum(abs(data.cons.hyb),2))/1e6);
% IDN
  if strcmp(param.from.machine,'JET')&(param.nombre.idn==16)
   fprintf(fid,'IDN (#%d) : Pmax = ',param.nombre.idn);fprintf(fid,'%4.1g ', ...
           [max(sum(real(data.cons.idn(:,1:8)),2))/1e6,max(sum(real(data.cons.idn(:,2:16)),2))/1e6]);	
   fprintf(fid,'MW, Etot = %.3g MJ\n',trapz(data.gene.temps,sum(real(data.cons.idn),2))/1e6);
  else
   fprintf(fid,'IDN (#%d) : Pmax = ',param.nombre.idn);fprintf(fid,'%4.1g ',max(real(data.cons.idn))/1e6);	
   fprintf(fid,'MW, Etot = %.3g MJ\n',trapz(data.gene.temps,sum(real(data.cons.idn),2))/1e6);
  end

% 4 - Gaz :
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  nb = fprintf(fid,'Gaz (Z,A) : ');
  lg = {'Maj','Min1','Min2','Imp1','Imp2'};
  for k =1:length(param.compo.z)
	if k <= length(lg)
		leg = lg{k};
	else
		leg = sprintf('Gaz%d',k);
	end
	if (param.compo.z(k) == 1) & (param.compo.a(k) ==1)
		nb = fprintf(fid,'%s = H, ',leg) + nb;
	elseif (param.compo.z(k) == 1) & (param.compo.a(k) ==2)
		nb = fprintf(fid,'%s = D, ',leg) + nb;
	elseif (param.compo.z(k) == 1) & (param.compo.a(k) ==3)
		nb = fprintf(fid,'%s = T, ',leg) + nb;
	elseif (param.compo.z(k) == 2) & (param.compo.a(k) ==3)
		nb = fprintf(fid,'%s = He3, ',leg) + nb;
	elseif (param.compo.z(k) == 2) & (param.compo.a(k) ==4)
		nb = fprintf(fid,'%s = He, ',leg) + nb;
	elseif (param.compo.z(k) == 6) & (param.compo.a(k) ==12)
		nb = fprintf(fid,'%s = C, ',leg) + nb;
	elseif (param.compo.z(k) == 8) & (param.compo.a(k) ==16)
		nb = fprintf(fid,'%s = O, ',leg) + nb;
	else
		nb = fprintf(fid,'%s = (%d,%d), ',leg,param.compo.z(k),param.compo.a(k)) + nb;
	end
	
	if (nb > 70) & (k <length(param.compo.z))
		fprintf('\n');
		nb =0;
	end
  end
%fprintf(fid,'Gaz :    Maj        Min1        Min2        Imp1        Imp2\n');
%fprintf(fid,'  Z :');fprintf(fid,'    %2d      ',param.compo.z);
%fprintf(fid,'\n  A :');fprintf(fid,'    %2d      ',param.compo.a);
  fprintf(fid,'\n<N0>:');
  for k =1:size(data.impur.impur,3)
   ss = sprintf('%8.3g',mean(data.impur.impur(:,1,k)));
   fprintf(fid,' %9.9s ',ss);
  end
  fprintf(fid,' m^-3');

  fprintf(fid,'\nZeff : <Zeffm>= %.3g, min(Zeffm)= %.3g, max(Zeffm)= %.3g\n', ...
        mean(data.cons.zeffm),min(data.cons.zeffm),max(data.cons.zeffm));
  fprintf(fid,'nH/nD: <nhnd> = %.3g, min(nhnd) = %.3g, max(nhnd) = %.3g\n', ...
        mean(data.cons.nhnd),min(data.cons.nhnd),max(data.cons.nhnd));
  fprintf(fid,'Nebar: <Nebar> = %.3g, min(Nebar) = %.3g, max(Nebar) = %.3g\t(1e19 m^-3)\n', ...
        mean(data.gene.nbar/1e19),min(data.gene.nbar/1e19),max(data.gene.nbar/1e19));

% 5- generation de courant :
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  fprintf(fid,'max de:   Ip       Icd       Ini       Iboot       Ihyb       Iidn       Ifci\n');
  fprintf(fid,'(en MA) %6.3g   %6.3g    %6.3g    %6.3g      %6.3g     %6.3g     %6.3g\n', ...
             max(data.gene.ip)/1e6,max(data.gene.icd)/1e6,max(data.gene.ini)/1e6,max(data.gene.iboot)/1e6, ...
	     max(data.gene.ihyb)/1e6,max(data.gene.iidn)/1e6,NaN);

% 5- configuration calcul :
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  fprintf(fid,'Mode   : Psi = %s, Pe = %s, Pion = %s, Ne = %s \n',zmode(data.mode.psi,langue),zmode(data.mode.pe,langue), ...
        zmode(data.mode.pion,langue),zmode(data.mode.nel,langue));
  spsi   = zlimite(data.mode.cons.psi,{'Ip','Vloop','Psi1','?'});
  spe    = zlimite(data.mode.cons.pe,{'Te1','qe1','Pe1','?'});
  spion  = zlimite(data.mode.cons.pion,{'Ti1','qion1','Pion1','?'});
  sne    = zlimite(data.mode.cons.ne,{'Ne1','ge1','?','?'});
  fprintf(fid,'Limite : Psi = %s, Pe = %s, Pion = %s, Ne = %s \n',spsi,spe,spion,sne);

  switch param.gene.cn
  case 0
    scn = 'implicite';
  case 0.5
    scn = 'C-N';
  case 1
    scn = 'explicite';
  otherwise
    scn = 'batard';
  end
  if param.gene.adiabatic  == 1
    sadia ='adiabatic';
  else
    sadia ='explicite';
  end
  fprintf(fid,'Solver : %s, amorti = %g, init = %s, d(adia) = %g, nmax = %d\n', ...
        scn,param.gene.amorti,sadia,param.gene.delta_adia,param.gene.nmax);

  switch param.gene.force
  case 0
    sfo = 'standard';
  otherwise
    sfo = 'forcee';
  end
  if param.gene.nonneg  == 0
    spos ='libre';
  else
    spos ='positif';
  end
  if param.gene.guido  == 0
    gg ='libre';
  else
    gg ='> 0';
  end
  fprintf(fid,'Special : convergence %s, sign(Pe, Pion & Ne) -> %s, J(0) %s \n', sfo,spos,gg);

  if param.gene.modecoef == 1
    scoef = 'convectif';
  else
    scoef = 'complet';
  end
  if param.gene.self == 1
    sself = 'selfconsitant';
  else
    sself = 'standard';
  end
  if param.gene.fast == 0
     sfast = 'selfconstitant';
  elseif param.gene.fast == 1
     sfast = 'rapide';
  else
    sfast = 'Jdiff';
  end	 	 
  fprintf(fid,'Eq : lambda= %g, #eq= %d, coef= %s, neo= %s, sources = %s\n',param.gene.lambda, ...
        param.gene.nbeq_mode,scoef,sfast,sself);
  switch param.gene.psiequi
  case 0
   sconv ='jamais';
  case 1
   sconv ='toujours';
  case 2
   sconv ='nonconvergence';
  otherwise
   sconv ='?';
  end
  fprintf(fid,'Equi : d(jmoy) = %g, amorti = %g, nmax = %d, equi -> diff = %s\n',param.gene.djmoy, ...
        param.gene.mjmoy,param.gene.nequi,sconv);

% manque sourcebord

  fprintf(fid,'Init : d(psi) = %g, nmax = %d\n',param.gene.dpsi_ini,param.gene.nequi_ini);
  fprintf(fid,'Critere : Psi = %g, Pe = %g, Pion = %g, Ne =%g\n',param.gene.critere.psi, ...
         param.gene.critere.pe,param.gene.critere.pion,param.gene.critere.ne);


% les modules
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  nom       = fieldnames(param.fonction);
  for k = 1:size(nom,1)
	if strcmp(nom{k},'mhd')
  
	else
		% paramatre du module
		parametre = getfield(param.cons,nom{k});
		% nom du module
		fonction  = getfield(param.fonction,nom{k});
		% mode associe s'il existe
		if isfield(data.mode,nom{k})
			mode   = getfield(data.mode,nom{k});
			switch nom{k}
			case {'glacon','plot'}
				smode = zlimite(mode,{'arret','marche','?','?'});
			otherwise
				smode  = zmode(mode,langue);
			end
		else
			smode   = '';
		end  
		
		% info complementaires  
		switch nom{k}
		case 'glacon'
			icomp = sprintf(', #%d',sum(data.cons.glacon(:)));
		case 'hyb'
			if isfield(parametre,'xdur')
				if getfield(parametre,'xdur') == 1
					icomp = ', xdur';
				else
					icomp = '';
				end
			end
		otherwise
			icomp = '';
		end
		
		% nomnbre de coupleur, injecteur ...
		if isfield(param.nombre,nom{k})
			nb     = getfield(param.nombre,nom{k});
		else
			nb     = 1;
		end
		if isempty(fonction)
			sparam = '';
		elseif ~isempty(parametre) 
			% parametre par defaut
			try
				info      = feval(fonction,nb);
				defaut    = info.valeur;
			catch
			   defaut    = [];    
			end
			% comparaison
			sparam = 'standard';
			pnom   = fieldnames(parametre);
			for l  = 1:size(pnom,1)
				pc  = getfield(parametre,pnom{l});
				if isfield(defaut,pnom{l})
					pd  = getfield(defaut,pnom{l});
					if length(pc) < length(pd)
					  pc(length(pd))=0;
					end
				elseif ischar(pc)
					pd = '';
				else
					pd =[];
				end
				
				if ischar(pc)
					if ~strcmp(pc,pd)
						sparam = 'personnalise';
						if length(icomp) < 16
							icomp = [icomp,', ',pnom{l}];
						elseif ~findstr(icomp,'...')
							icomp = [icomp,', ...'];
						end
					end
				elseif  ~isempty(pc) & ~isempty(pd) & isnumeric(pc)
					if ~all(size(pc) == size(pd))
						sparam = 'personnalise';
						if length(icomp) < 16
							icomp = [icomp,', ',pnom{l}];
						elseif ~findstr(icomp,'...')
							icomp = [icomp,', ...'];
						end
					elseif ~all(all(pc==pd))
						sparam = 'personnalise';
						if length(icomp) < 16
							icomp = [icomp,', ',pnom{l}];
						elseif ~findstr(icomp,'...')
							icomp = [icomp,', ...'];
						end
					end
				end
			end 
		else
			sparam = 'sans';
		end
		
		if isempty(fonction)
			% rien
		elseif isempty(smode)
			fprintf(fid,'%8.8s : %12.12s (%12.12s)%s\n',nom{k},fonction,sparam,icomp); 
		else
			fprintf(fid,'%8.8s : %12.12s (%12.12s) -> %8.8s %s\n',nom{k},fonction,sparam,smode,icomp); 
		end
	end	 
  end

% coef de transport
  fprintf(fid,'-------------------------------------------------------------------------------\n');
% a - les multiplicateur
  fprintf(fid,'neo *    : ');
  neonom = fieldnames(param.cons.neomulti);
  for k =1:length(neonom)
	if k <length(neonom)
		fprintf(fid,'%2s = %4.2g, ',neonom{k},getfield(param.cons.neomulti,neonom{k}));
	else
		fprintf(fid,'%2s = %4.2g',neonom{k},getfield(param.cons.neomulti,neonom{k}));
	end	
  end
  fprintf(fid,'\n');
  coefnom = fieldnames(data.coef);
  smode ={};
  for k =1:length(coefnom)
	if isfield(data.mode,coefnom{k})
		mode   = getfield(data.mode,coefnom{k});
		smode{k}  = zmode(mode,langue);
		
	else
		smode{k}='';
	end
  end
  if ~isempty(strmatch('zero',smode))
	nb = fprintf(fid,'a zero   : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'zero')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('donnees',smode))
	nb = fprintf(fid,'donnees  : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'donnees')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('calcule',smode))
	nb = fprintf(fid,'calcule  : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'calcule')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('complexe',smode))
	nb = fprintf(fid,'complexe : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'complexe')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
% les modes
  fprintf(fid,'-------------------------------------------------------------------------------\n');
  coefnom = fieldnames(data.mode);
  smode ={};
  exclusion = {'pe','pion','nel','psi','premiertemps'};
  for k =1:length(coefnom)
	if ~isstruct(getfield(data.mode,coefnom{k}))
		if ~isfield(param.fonction,coefnom{k}) &  ...
		   ~isfield(data.coef,coefnom{k}) &  ...	
		   isempty(strmatch(coefnom{k},exclusion,'exact'))
			mode   = getfield(data.mode,coefnom{k});
			smode{k}  = zmode(mode,langue);
		else
			smode{k}='';
		end
	else
			smode{k}='';
	end
  end
  if ~isempty(strmatch('zero',smode))
	nb = fprintf(fid,'a zero   : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'zero')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('donnees',smode))
	nb = fprintf(fid,'donnees  : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'donnees')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('calcule',smode))
	nb = fprintf(fid,'calcule  : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'calcule')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
  if ~isempty(strmatch('complexe',smode))
	nb = fprintf(fid,'complexe : ');
	for k = 1:length(coefnom)
		if strcmp(smode{k},'complexe')
			nb = fprintf(fid,'%s ',coefnom{k});
		end
		if nb > 75
			nb = fprintf(fid,'\n           ')-1;
		end
	end
	fprintf(fid,'\n');
  end
% option de creation
  fprintf(fid,'-------------------------------------------------------------------------------\n');
% parametre par defaut
  if ~isfield(param.from,'createur')
    param.from.createur ='';
  end
    
  if ~isempty(param.from.createur)
	createur  = param.from.createur;
	if strcmp(createur,'zvplusacces')
	  createur = 'zplusacces';
	end
	info      = feval(createur);
	defaut    = info.valeur;
	parametre = param.from.option;
	nb = fprintf(fid,'createur : %s -> ',createur);
	% comparaison
        if ~isempty(parametre)
	   pnom   = fieldnames(parametre);
	   for l  = 1:size(pnom,1)
		   pc  = getfield(parametre,pnom{l});
		   if isfield(defaut,pnom{l})
			   pd  = getfield(defaut,pnom{l});
		   elseif ischar(pc)
			   pd = '';
		   else
			   pd =[];
		   end
		   if ischar(pc)
			   if ~strcmp(pc,pd)
				   nb = nb + fprintf(fid,'%s = %s, ',pnom{l},pc);
			   end
		   elseif isnumeric(pc)
			   if  ~isempty(pc) & ~isempty(pd) & ~all(all(pc==pd))
				   nb = nb + fprintf(fid,'%s = %g, ',pnom{l},pc);
			   end
		   end
		   if nb > 70
			   fprintf(fid,'\n');
			   nb = 0;
		   end
	   end	
        end            
  else
	nb = fprintf(fid,'createur : ???');
  end 
             
             
% sortie 
%x -fichier
  fprintf(fid,'\n-------------------------------------------------------------------------------\n');
  fprintf(fid,'fichier associe :\t%s\n',fichier);
  [s,t]=unix(['ls ',param.gene.origine,'.*']);
  if s ~= 0
	fprintf(fid,'le fichier, definit dans param.gene.origine, n''existe plus !\n');
  elseif strcmp(param.gene.filetype,'resultat')
	fprintf(fid,'fichier de sortie : \n  %s\n',param.gene.origine);
  else
	fprintf(fid,'fichier d''entree : \n  %s\n',param.gene.origine);
  end
  [s,t]=unix(['ls ',param.gene.file,'.*']);
  if s ~= 0
    fprintf(fid,'le fichier, definit dans param.gene.file, n''existe plus !\n');
  else
    fprintf(fid,'fichier de sortie : \n  %s\n',param.gene.file); 
  end

% recherche du status
  if isempty(param.gene.rapsauve)
	fprintf(fid,'pas de sauvegarde rapide definie\n');
  else
	rap = param.gene.rapsauve;
	ind = max(findstr(rap,'data/'));
	if ~isempty(ind)
		rap = strcat('...',rap((ind + 4):end));
	end
	
	[s,t] = unix(['ls -1 ',param.gene.rapsauve,'*']);
	if s ~= 0
		fprintf(fid,'nom pour la sauvegarde rapide (vide) : \n  %s\n',rap);
	else
		nb = size(tseparec(t),1);
		fprintf(fid,'nom pour la sauvegarde rapide (%d temps a reconstruire): \n  %s\n',nb,rap);
	end
  end
% fermeture du fichier si pas sortie standart
  if fid > 2
    fclose(fid);
  end
% affichage du resume si necessaire
  if nargin  == 4
    [s,message] = unix(['cat ',file]);
    if s ~= 0
      disp('Probleme de lecture du resume :')
		disp(message);
		message ='';
    end
    [s,t] = unix(['rm -f ',file]);
    if s ~= 0
        disp('Probleme pour effacer le fichier temporaire du resume :')
		disp(t)
    end
  elseif edit == 1
     if efface == 1
       [s,t] = unix([getappdata(0,'editeur'),' ',file,' &']);
       if strcmp(computer,'ALPHA')
         [s,t] = unix(['sleep 60;rm -f ',file,'&']);
       end
     else
       [s,t] = unix([getappdata(0,'editeur'),' ',file,' &']);
     end
     if s ~= 0
       disp('Probleme unix :')
       disp(t)
       cr = 107;
       return
     end 
  elseif isempty(getenv('DISPLAY')) & ~isempty( param.from.creation.user)
  % envoi par mail si batch
    [s,message] = unix(['cat ',file]);
    cr_mail=zmail(strtok(param.from.creation.user),sprintf('Resume Cronos : %s',file),message);
  end
end
% format affichage date
function s = zdate(vec)

try
    n = datenum(vec(1),vec(2),vec(3),vec(4),vec(5),vec(6));
    s = datestr(n);
catch
    s ='???';
end

% test des modes
function s =zmode(mode,langue)

if strcmp(langue,'anglais')
  if all(mode == 0)
   s = 'null';
  elseif all(mode == 1)
   s = 'input';
  elseif all(mode == 2)
   s = 'calculated';
  elseif all(mode == 3)
   s = 'recopy';
  else
   s = 'complex';
  end
else
  if all(mode == 0)
   s = 'zero';
  elseif all(mode == 1)
   s = 'donnees';
  elseif all(mode == 2)
   s = 'calcule';
  elseif all(mode == 3)
   s = 'recopie';
  else
   s = 'complexe';
  end

end
% test des modes
function s =zmode2(mode,langue)

if all(mode == 0)
	s = 'null';
elseif all(mode == 1)
	s = 'input';
elseif all(mode == 2)
	s = 'calculated';
elseif all(mode == 3)
	s = 'density total controled';
elseif all(mode == 4)
	s = 'density total controled + shape law in ngr/nbar ';
elseif all(mode == 5)
	s = 'density total controled + shape Wiesen';
else
	s = 'complex';
end


% test des conditions aux limites
function s =zlimite(mode,info)

if all(mode == 0)
   s = info{1};
elseif all(mode == 1)
   s = info{2};
elseif all(mode == 2)
   s = info{3};
elseif all(mode == 3)
   s = info{4};
else
   s = 'complexe';
end
