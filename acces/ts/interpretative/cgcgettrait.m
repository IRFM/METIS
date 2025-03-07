% CGCGETTRAIT rammene dans une structure ou dans l'envirronement d'appel les donnees d'un traitement
%---------------------------------------------------------------------------------------------------
% fichier cgcgettrait.m -> cgcgettrait, zzdia
%
% Cette fonction ramene dans l'environnement de
% l'appelant ou dans une structure, les donnees
% produite par un traitement (ou un diagnostic).
% 
% Si le nom est precede ou suivi d'un ";"
% interprete cela comme une consigne ou 
% un asservissement
%
% Le nom des variables est celui donne dans la base 
% moins le prefixe ('s','g','r','sc'). 
% Pour les traitement suivant, la racine du nom du 
% traitement est aussi retiree :
% 
%  tprof      -> prof
%  tcoupol    -> coup
%  tboostrap  -> boot
%  tcronos    -> cron
%  tbragg     -> brag
%  tfce       -> fce
%
% valeurs speciales :
%     
%   nomtrait = 'loco' ramenne les donnees du fichier loco
%
% syntaxe : 
%
%    1- rerour dans une structure : 
%    
%            [data,sout] = cgcgettrait(numchoc,nomtrait,{subs});
%    
%    2- retour dans l'environnement de l'appelant :
%
%            cgcgettrait(numchoc,nomtrait,{subs});
%            
%            
% entrees :
% 
%     numchoc  = numero du choc.numero occurrence
%     nomtrait = nom du traitement (ou du diagnostic)
%     subs     = racine de nom de traitement a retirer
%                (par defaut donne par le nom du traitement)
%
% sortie : 
% 
%     data  = structure de donnees de sortie simple 
%             (contient que les donnees + une variable
%             de temps 'times' et une donnee d'espace
%             rhofit). C'est les champs de cette structure 
%             qui sont retourner dans l'envirronnement de 
%             l'appelant. 
%            
%     sout  = structure de sortie des donnees, ou chaque donnee est 
%             associe a une sous structure de format :
%                    nomdonnee.data           -> les donnees
%                    nomdonnee.temps          -> la coordonnee de temps 
%                    nomdonnee.espace         -> la coordonnee d'espace
%                    nomdonnee.information    -> la struture d'information
%                    nomdonnee.commentaire    -> le commentaire associe a la donnee
%               
%  remarque :
%    
%    1- la sous structure info :    (de tsbase_cert)
%    
%         information.cert    =  certification de la donnee
%         information.vers    =  version base du traitement
%         information.date    =  date d'ecriture dans la base
%         information.heure   =  heure d'ecriture dans la base
%         information.unix    =  unite du temps
%         information.uniy    =  unite des donnees
%         information.uniz    =  unite d'espace
%         
%    2- le signal commentaire est copie dans le champ texte_complet. De plus,
%       pour les supertraitement il est decode dans une sous structure 'info'.
%       
%    3- si nomtrait se termine par '@', la variable times de l'appelant devient la
%       variable 'times'  de la structure data. 
%       
%    4- si les temps des donnees sont differents, les donnees de la structure data sont
%       reechantillonees sur 'times' des qu'elle est definie
%
% fonction ecrite par J-F Artaud, poste 46-78
% version 1.0, derniere modification le 25/10/1999
% modification :
% 
% 
% 
%---------------------------------------------------------------------------------------------------
%
function [data,sout] = cgcgettrait(numchoc,nomtrait,subs)

% sorties a vide
sout =[];
data.times=[];
data.rhofit=[];

% test des entrees
if nargin <2
	error('Il faut donnee le numero de choc et le nom du traitement');
elseif isempty(numchoc)
	error('numero de choc vide');
elseif isempty(nomtrait)
	error('nom du traitement vide');
end

if nargin <3
   subs='';
end
 
% si reechantillonage
if iscell(nomtrait)

    nom1 = nomtrait{1};
	 if nom1(end)=='@'
    	nomtrait{1}=nom1(1:(end-1));
    	recupere =1;
	else
   	 recupere =0;
	end
   nomtrait = strvcat(nomtrait);
else
	if nomtrait(end)=='@'
   	 nomtrait=nomtrait(1:(end-1));
    	recupere =1;
	else
   	 recupere =0;
	end
end
% detection des cas particulier
particulier = 'loco';
if ~isempty(strmatch(deblank(lower(nomtrait)),particulier,'exact'))
   
   % acces au donnees depuis un fichier
   
   if strcmp(deblank(lower(nomtrait)),'loco')
        
        nomfile = strcat('/usr/drfc/cgc/st/data/',num2str(numchoc),'/loc.mat');
        
        % lecture des donnees
        try
            lo = load(nomfile);
        catch
            warning('Pas de donnees pour ce choc')
            return
        end
        
        % la base temps
        if recupere == 1
           try
	      data.times=evalin('caller','times');
           catch
	      data.times =lo.loc.TIME;
           end
        else
	   data.times =lo.loc.TIME;
        end
        data.rhofit = lo.loc.RHO';
        t =lo.loc.TIME;
        v = lo.loc.RHO';
        
        % eclatement de la structure
        donnees = struct2cell(lo.loc);
        nom     = fieldnames(lo.loc);
   
        % boucle de contstruction des  srcutures sout et data
        for k =1:length(nom)
            d= donnees{k};
%            sout = setfield(sout,strcat(nom{k},'.commentaire'),'non demande ou vide');
%            sout = setfield(sout,strcat(nom{k},'.information'),'');
%            sout = setfield(sout,strcat(nom{k},'.data'),d);
            sout = setfield(sout,nom{k},'commentaire','non demande ou vide');
            sout = setfield(sout,nom{k},'information','');
            sout = setfield(sout,nom{k},'data',d);
            if size(d,2) > 1 
%                sout = setfield(sout,strcat(nom{k},'.espace'),v);
                sout = setfield(sout,nom{k},'espace',v);
            else
%                sout = setfield(sout,strcat(nom{k},'.espace'),[]);
                sout = setfield(sout,nom{k},'espace',[]);
            end
            if length(t) == size(d,1)
%               sout = setfield(sout,strcat(nom{k},'.temps'),t);
               sout = setfield(sout,nom{k},'temps',t);
               if ~isequal(t,data.times) & (size(t,1) > 2)
           	   try
           	      d = tsample(d,t,data.times,'fen');
                      d = zinter(d,t,data.times);
           	   catch 
           	      warning('Pb tsample');
           	   end
               end
            elseif length(lo.loc.TTIME) == size(d,1)
%               sout = setfield(sout,strcat(nom{k},'.temps'),lo.loc.TTIME);
               sout = setfield(sout,nom{k},'temps',lo.loc.TTIME);
               try
           	  d = tsample(d,lo.loc.TTIME,data.times,'fen');
                  d = zinter(d,lo.loc.TTIME,data.times);
               catch 
           	  warning('Pb tsample');
               end
            else
%               sout = setfield(sout,strcat(nom{k},'.temps'),[]);
               sout = setfield(sout,nom{k},'temps',[]);
               d=d(:)';
            end
            data = setfield(data,nom{k},d);
        end
        
        % donnees complementaires
        sout.job = lo.job;
        sout.don = lo.don;
        sout.info.num = lo.num;
        sout.info.clock = lo.clock_info;
        sout.info.occurrence = lo.occurrence;
        sout.info.user = lo.user;
        sout.info.system = lo.system;
        data.job = lo.job;
        data.don = lo.don;
        data.info.num = lo.num;
        data.info.clock = lo.clock_info;
        data.info.occurrence = lo.occurrence;
        data.info.user = lo.user;
        data.info.system = lo.system;
   end
else

    nomtrait =lower(deblank(nomtrait));
    modedia  = 0;
    if  all(size(nomtrait) >1)
	    % liste de signaux
	    liste = nomtrait;
	    commentaire ='';
    elseif ~isempty(findstr(nomtrait,';'))
            % cas des consignes
 	    nomtrait = strrep(nomtrait,';','');
            [liste,commentaire] = zzcons(fix(numchoc),nomtrait);
    elseif nomtrait(1)=='d'
	    % cas du diagnostic
	    numchoc=fix(numchoc);
	    [liste,commentaire] = zzdia(numchoc,nomtrait);
	    modedia = 1;
    elseif ((tsexist(numchoc,nomtrait) == 0) &(nomtrait(1)~='t'))| ...
           ((nomtrait(1)=='s')|(nomtrait(1)=='g')|(nomtrait(1)=='r'))
	    % cas du signal simple
	    liste = nomtrait;
	    commentaire ='';
    else
       [liste,commentaire] = tsbase(nomtrait);
    end

    if isempty(liste)
	    return
    end

    % racine dans le nom
    if ~isempty(subs)
         % on ne fait rien
    elseif strcmp(nomtrait,'tprof')
          subs  = 'prof';
    elseif strcmp(nomtrait,'tcoupol')
          subs  = 'coup';
    elseif strcmp(nomtrait,'tbootstrap')
          subs  = 'boot';
    elseif strcmp(nomtrait,'tcronos')
          subs  = 'cron';
    elseif strcmp(nomtrait,'tbragg')
          subs  = 'brag';
    elseif strcmp(nomtrait,'tfce')
          subs  = 'fce';
    else
          subs ='';
    end

    if recupere == 1
        try
	    data.times=evalin('caller','times');
        catch
	    data.times =[];
        end
    end


    % boucle sur les signaux
    for k = 1:size(liste,1)
	    d=[];
	    t=[];
	    v=[];
	    c='';
	    decode =0;
	    signal = lower( deblank(liste(k,:)));
            if ~isempty(findstr(signal,';'))
                    % cas des consignes
                    aaa = tsmat(fix(numchoc),strcat(deblank(liste(k,:)),';VALEURS'));
                    ind = findstr(signal,';');
                    signal = signal(ind(1):end);
                    signal = strrep(signal,'+','plus');
                    signal = strrep(signal,'&','and');
                    signal = strrep(signal,'-','moins');
                    signal = strrep(signal,'/','slash');
                    signal = strrep(signal,'*','star');
                    if ~isempty(aaa)
                         d = aaa(:,2);
                         t = aaa(:,1);
                    end
                    if isempty(data.times)
                        % creation de la base temps consignes
                        [void,times] = tsbase(fix(numchoc),'simag');
                         data.times = times;
                    end
	    elseif lower(signal(1)) == 'g'
		    [d,t,v,c]=tsbase(numchoc,signal);
	    elseif strcmp(lower(signal(1:2)),'sc') 
		    if isempty(strmatch(lower(signal),'scionbrag'))
		    		[d,t,c,vv]=tsbase(numchoc,signal);	
			 else
					d=[];t=[];c=[];v =[];
					disp('Probleme base de donnees : scionbrag illisible')
			 end
		 else
		    [d,t,c]=tsbase(numchoc,signal);		
	    end	
	    
            if signal(1) == 'g'
                signal =signal(2:end);
            elseif ischar(d)  & strcmp(signal(1:2),'sc')
					 if strcmp(lower(signal),'scticxs')
					  	signal =signal(2:end);
					 elseif ~strcmp(subs,'cron') 
                				signal =signal(3:end);
				         
					 else
                				signal =signal(2:end);
					 end
                t=[];
                if ~isempty(subs)
					    if ~strcmp(subs,'cron') & ~strcmp(subs,'brag') & ~strcmp(subs,'fce') 
                   	decode =1;
						 end
                end
            elseif  signal(1) == 'r'
                signal =signal(2:end);
                t=[];
            else
                signal =signal(2:end);
            end

	    signal=strrep(signal,subs,'');

	    if decode == 0
		         % si signal vide
					if isempty(signal) & strcmp(subs,'cron')
					  	 signal = 'resume'
				   elseif isempty(signal)
					  	 signal = 'commentaire';
					end
               % remplissage de la structure de sortie
%               sout = setfield(sout,strcat(signal,'.data'),d);
%               sout = setfield(sout,strcat(signal,'.temps'),t);
%               sout = setfield(sout,strcat(signal,'.espace'),v);
               sout = setfield(sout,signal,'data',d);
               sout = setfield(sout,signal,'temps',t);
               sout = setfield(sout,signal,'espace',v);

               % certification
               warning off
               [info.cert,info.vers,info.date,info.heure, ...
                info.unix,info.uniy,info.uniz,info.nomdon, ...
                                      info.type]=tsbase_cert(c);
                info.num = numchoc;
               warning on
%               sout = setfield(sout,strcat(signal,'.information'),info);
               sout = setfield(sout,signal,'information',info);

               % commentaire
               if ~isempty(commentaire)
%                   sout = setfield(sout,strcat(signal,'.commentaire'),commentaire(k,:));
                   sout = setfield(sout,signal,'commentaire',commentaire(k,:));
               else
%                   sout = setfield(sout,strcat(signal,'.commentaire'),'non demande ou vide');
                   sout = setfield(sout,signal,'commentaire','non demande ou vide');
               end

               % struture data 
	       if modedia  == 1
	           % correction des basestemps
	           ddt = diff(t,1,1);
		   ddtl = ddt(:);
		   dtm = mean(ddtl(ddtl>0));
		   ind = find(ddt <= 0);
		   if ~isempty(ind)
		     t(ind - 1) = t(ind) +dtm;
		   end
	       end
               if isempty(data.times) & ~isempty(t)
           	    data.times = t(:,1);
               elseif ~isempty(t)
           	    if ~isequal(t,data.times) & (size(t,1) > 2)
           		    try
			       if length(data.times) < 3
				   d = interp1(d,t,data.times,'nearest','extrap');
                               elseif length(data.times) > (size(t,1)*3)
           		           d = tsample(d,t,data.times,'fel');
                               else
           		           d = tsample(d,t,data.times,'fen');
                               end
                               d = zinter(d,t,data.times);
           	            catch 
           	               warning('Pb tsample');
           	            end
           	    end
               end
               if isempty(data.rhofit) & ~isempty(v)
                  data.rhofit = v;
               elseif ~isempty(v)
           	    if ~isequal(v,data.rhofit)
%           		    data = setfield(data,strcat('r.',signal),v);
           		    data = setfield(data,'r',signal,v);
           	    end
               end
               data = setfield(data,signal,d);
            else
               data = setfield(data,'texte_complet',d);
               sout = setfield(sout,'texte_complet',d);
   	       d=tseparec(d);
     	       for l =1:size(d,1);
   		    			[nom,ch]=strtok(d(l,:),' ');
 	            		signal = strcat('info.',nom);
							if isempty(ch)
							    ch = 'vide';
								 signal(end)=[];   % pour enlever le point a la fin (FI, 25/07/02)
							end
%                    data = setfield(data,signal,ch);
%                    sout = setfield(sout,signal,ch);
						  eval(sprintf('data.%s=ch;',signal));
						  eval(sprintf('sout.%s=ch;',signal));
               end

            end


    end
end
    
if nargout == 0 & ~isempty(data)
    donnees = struct2cell(data);
    nom     = fieldnames(data);
    for k = 1:length(nom)
       assignin('caller',nom{k},donnees{k});
    end
    clear data sout
end


% fonction zzdia d'extraction d'es informations sur un diagnostic
function [l,c]=zzdia(n,p)

l=[];
c=[];
p =upper(deblank(p));
lt=tsmat(fix(n));
if isempty(lt)
	return
end

ind1= min(strmatch(p,lt));
if isempty(ind1)
	return
end
indd= strmatch('D',lt(:,1));
ind2=min(indd(indd>ind1));
if isempty(ind2)
        indd= strmatch('EXP',lt);
        if indd(1) > ind1
            ind2 = indd(1);
        else
            ind2 = length(lt);
        end
end
v=(ind1+1):(ind2-1);
l=lt(v,:);
c=lt(v,:);


% met des NaN hors de l'intervalle de validite
function d = zinter(d,t,te);

      ind =  find((te < min(t(:)))|(te > max(t(:))));
      d(ind,:) = NaN .* ones(length(ind),size(d,2));
        
% retourne la liste des consignes
function  [liste,commentaire] = zzcons(numchoc,nomtrait)

liste        = '';
commentaire  = '';
% lecture de la liste des asservissement
%requete ='select nom_objet, producteur from objet where type_objet=21 and ensemble = ''choc'' order by producteur';
requete ='select distinct nom_objet, producteur from objet where type_objet=21  order by producteur';
[cr,temps,requete,texte,txt,nb_para,nom_para,cons,prod]=matlsql('','',requete,ones(40,1));

% recherche d'un producteur
ind = strmatch(lower(nomtrait),lower(prod),'exact');
if ~isempty(ind)
     for k =1:length(ind)
          nomc  = strcat(prod(ind(k),:),';',cons(ind(k),:));
          if isempty(strmatch(nomc,liste,'exact'))
               liste = strvcat(liste,nomc);
          end
     end
else
    ind = min(strmatch(lower(nomtrait),lower(cons),'exact'));
    if isempty(ind)
        ind = strmatch(lower(nomtrait),lower(cons));
    end
    if isempty(ind)
        return
    end
     for k =1:length(ind)
          nomc  = strcat(prod(ind(k),:),';',cons(ind(k),:));
          if isempty(strmatch(nomc,liste,'exact'))
               liste = strvcat(liste,nomc);
          end
     end
end
commentaire  = liste;

