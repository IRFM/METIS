% MATLSQL	execute une requete sql sur une base de donnees
%-----------------------------------------------------------------
%
%	fonction matlsql.m
%
%	cette fonction execute depuis matlab une requete sql en
%	direction du serveur Ingres choisi sur la base de son choix.
%
%	la requete est un texte sql sans le ';' a la fin  
%	(il est rajouter afin que l'execution soit garantie)
%
%	elle remplace la fonction tssql.m (plus rapide et plus complete sur le decodage)
%
%	syntaxe :
%		
%		[cr,temps,requete,texte,txt,nb_para,nom_para,s1, ... ,{s40}]=tssql(serveur,base,requete,formats);
%
%	entree :
%
%		serveur		=	nom du serveur Ingres (defaut = 'rigel')
%
%		base		=	nom de la base de donnees (defaut = 'arcad') 
%
%		requete		=	vecteur ou matrice de texte matlab
%						contenant la requete
%		format_s		=	indique le formats des sorties s1 ...s40.
%						cette variable est optionnelle. Par defaut le format
%						est detecter automatiquement (ca peut etre lent).
%						format_s doit avoir une longueur de 40.
%						format_s(k) =	0	-> le format de sk est chosis automatiquement
%										1	-> le format de sk est alphanumerique
%
%	sortie :
%
%		cr	=	compte rendu : 0-> ok 
%
%		temps	= temps mis pour executer la requete
%
%		requete	= texte complete de la requete transmise aui moniteur sql
%		
%		texte	= reponse de la base 
%
%		txt		= partie utile de la reponse
%
%		nb_para = nombre de parametres extraits de la base
%
%		nom_para = nom des parametres 
%
%		s1 ... s40	=	variables extraites de la reponse du moniteur SQL
%
%	fonction ecrite par J-F Artaud, poste 46-78
%	version 1, derniere modification le 24/03/98
%-----------------------------------------------------------------
%
function [	cr,temps,requete,texte,txt,nb_para,nom_para, ...
			s1,s2,s3,s4,s5,s6,s7,s8,s9,s10, ...
			s11,s12,s13,s14,s15,s16,s17,s18,s19,s20, ...
			s21,s22,s23,s24,s25,s26,s27,s28,s29,s30, ...
			s31,s32,s33,s34,s35,s36,s37,s38,s39,s40 ...
							]=matlsql(serveur,base,requete,format_s)

	%
	% les sorties
	%
	cr=-1;
	texte='';
	nb_para=0;
	nom_para='';
	s1=[];s2=[];s3=[];s4=[];s5=[];s6=[];s7=[];s8=[];s9=[];s10=[];
	s11=[];s12=[];s13=[];s14=[];s15=[];s16=[];s17=[];s18=[];s19=[];s20=[];
	s21=[];s22=[];s23=[];s24=[];s25=[];s26=[];s27=[];s28=[];s29=[];s30=[];
	s31=[];s32=[];s33=[];s34=[];s35=[];s36=[];s37=[];s38=[];s39=[];s40=[];
	maxnbch=40;
	txt='';
	%
	% recupere les serveurs et base par defaut
	%
	nom_base=getenv('NOM_BASE');
	if ~isempty(nom_base),
		ind=find(nom_base==':');
		if length(ind)==2,
			serveur_def=nom_base(1:(ind(1)-1));
			base_def=nom_base((ind(2)+1):length(nom_base));
		else
			serveur_def='deneb';
			base_def='arcad';
		end
	else
		serveur_def='deneb';
		base_def='arcad';
	end
	%
	% gestion des arguments
	%
	if nargin <1,
			serveur=serveur_def;
	end
	if nargin <2,
		base=base_def;
	end
	if nargin <3,
		requete='help';
	end
	if nargin <4,
		format_s=zeros(maxnbch,1);
	end
	%
	% variables vides
	%
	if isempty(serveur),
		serveur=serveur_def;
	end
	%
	if isempty(base),
		base=base_def;
	end
	%
	if isempty(requete),
		requete='help';
	end
	%
	if length(format_s)~=40,
		format_s=zeros(maxnbch,1);
	end
	%
	% enleve les blancs initiles
	%
	serveur(serveur==' ')=[];
	base(base==' ')=[];
	%
	% requete sur une ligne
	%
	if size(requete,1)>1,
		texte=[];
		for k=1:size(requete,1),
			texte=[texte,requete(k,:),10];
		end
		requete=texte;
	end
	%
	% recherche des instruction "in"
	%
	resquete=[requete,'     '];
	ind_in=findstr(lower(requete),' in');
	if ~isempty(ind_in)
		ind_in=[ind_in,length(requete)];
		for k_in =1:(length(ind_in)-1)
			txt_in=requete((ind_in(k_in)+3):(ind_in(k_in+1)-1));
			if ~isempty(txt_in)
				ind_p1=min(find(txt_in =='('));
				ind_p2=min(find(txt_in ==')'));
				if (~isempty(ind_p1))&(~isempty(ind_p2))& ...
					((txt_in(1)==' ')|(txt_in(1)=='('))
					ind_vir=find(txt_in==',');
					if length(ind_vir)>30
						fprintf('La requete contient un operateur de predicat "in" suivi d''une liste contenant plus de 30 elements.\n');
						fprintf('du fait d''un probleme Ingres/OSF, cette requete ne peut pas etre executee.\n');
						fprintf('vous devez scinder la requete ou utiliser une table temporaire.\n');
						fprintf('contactez votre administrateur de base de donnees\n');
										[cr,t]=unix(['Mail -s "requete avec IN illegale" eymin <<@@',10, ...
													'le user ',getenv('USER'), ...
													' a essaye d''executer une requete utilisant un IN avec', ...
													int2str(length(ind_vir)),' elements dans la liste',10, ...
													'La requete debute par :',10, ...
													requete(1:min(76,length(requete))),' ...' , 10,'.',10,'@@',10 ]);
	
						return
					end
				end
			end
		end
	end
	%
	% composition de la commande
	%
	requete((requete <32)&(requete~=10))=[];
	requete=sprintf('<<@@\n%s;\n',requete);
	%
	% commande  sql pour l'execution et la fin
 	%
	requete=[requete,sprintf('%s\n','\p')];
	requete=[requete,sprintf('%s\n','\go')];
	requete=[requete,sprintf('%s\n@@\n','\quit')];
	%
	% execution
	%
	fprintf('%s','execution de la requete');	
	temps=clock;
	%[cr,texte]=unix(['cd /usr/local/ingres/bin/;', ...
	%				 'sql -c1 -t1 -f4n11.3 -f8n15.3 ', ...
	%				 serveur,'::',base, ...
	%				 ' ',requete]);
    % test where is tm program
    loc_ingres = getenv('II_SYSTEM');
    if ~isempty(loc_ingres)
        cmd = fullfile(loc_ingres,'ingres','bin','tm'); 
    else
        [cr,texte] = zunix('which tm');
        if cr ~= 0
            cmd = '/opt/Ingres/IngresI0/ingres/bin/tm';
        else
            cmd = texte(texte > 32);
            cmd = strtrim(texte);
        end
    end
    %sprintf('exec %s -qSQL -c1 -t1 -f4n11.3 -f8n15.3 %s::%s %s', ...
	%				 cmd,serveur,base,requete)
	[cr,texte]=zunix(sprintf('exec %s -qSQL -c1 -t1 -f4n11.3 -f8n15.3 %s::%s %s', ...
					 cmd,serveur,base,requete));
	temps=etime(clock,temps);
	%
	if cr~=0
		fprintf('\n%s\n',['erreur dans l''execution de la requete sur ', ...
							serveur,'::',base]);
		return
	elseif isempty(texte),
		fprintf('%s\n','Pas de donnees');
		return
	else
		fprintf('%s\n',', puis decodage');
	end
	%
	% extraction de la partie utile de la reponse
	%
	txt=texte;
	txt(txt==13)=[];			% EOF
	%
	%recherche de la reponse
	%
	ind1=min(findstr(txt,'* Executing'));
	if isempty(ind1),
		ind1=1;
	end
	ind2=max(findstr(txt,'row'));
	if isempty(ind2),
		ind2=max(findstr(txt,'continue'));
	end	
	if isempty(ind2),
		ind2=length(txt);
	end
	%
	% separes les lignes
	%
	txt=tseparec(txt(ind1:ind2),10);
	if size(txt,1) >=3,
		txt=txt(2:(size(txt,1)-1),:);
	else
		txt=[];
	end
	%
	% enleve l'encadrement
	%
	if ~isempty(txt)
		ind=((txt=='+')|(txt=='-')|(txt==' '));
		if any(any(ind))
			li=1:size(txt,1);
			li=li(~all(ind')');
			txt=txt(li,:);
		end
	end
	%
	% enleve les ligne vides
	%
	if ~isempty(txt)
		ind=(txt==' ');
		if any(any(ind))
			li=1:size(txt,1);
			li=li(~all(ind')');
			txt=txt(li,:);
		end
	end
	%
	% choix du mode de decodage
	%
	if isempty(txt),
		disp('Pas de donnees utiles');
	elseif all((txt(:,1)=='|')|(txt(:,1)==' ')),
		%
		% recherche des caracteres '\x'
		%
		ind=find(any(txt'=='\'));
		if ~isempty(ind),
			for k=ind,
				litxt=txt(k,:);
				nbrep=sum(litxt=='\');
				litxt=strrep(litxt,'\a',char(7));
				litxt=strrep(litxt,'\b',char(8));
				litxt=strrep(litxt,'\t',char(9));
				litxt=strrep(litxt,'\n',char(10));
				litxt=strrep(litxt,'\v',char(11));
				litxt=strrep(litxt,'\f',char(12));
				litxt=strrep(litxt,'\r',char(13));
				%
				litxt=strrep(litxt,'\A',char(7));
				litxt=strrep(litxt,'\B',char(8));
				litxt=strrep(litxt,'\T',char(9));
				litxt=strrep(litxt,'\N',char(10));
				litxt=strrep(litxt,'\V',char(11));
				litxt=strrep(litxt,'\F',char(12));
				litxt=strrep(litxt,'\R',char(13));
				%
				litxt=strrep(litxt,'\\','\');
				litxt=strrep(litxt,'\"','"');
				litxt=strrep(litxt,'\''','''');				
				litxt=strrep(litxt,'\0',char(0));
				%
				if any(litxt=='\'),
					litxt(litxt=='\')=[];
					disp('Warning: Le decodage est peut-etre mauvais, verifier le nombre de parametres en sortie')
				end
				comp=char(32*ones(1,nbrep));
				txt(k,:)=[litxt,comp];
			end
		end
		%
		% recuperation de donnees sur select
		% les indice des collonnes
		%
		ind=find(all(txt=='|'));
		nb_para=length(ind)-1;
		nbli=size(txt,1);
		%
		% boucle sur les champs
		%
		for k=1:min(nb_para,maxnbch),
			%
			% la variable au format texte
			%
			nom_para(k,:)=setstrl(txt(1,(ind(k)+1):(ind(k+1)-1)),32);
			champ=txt(2:nbli,(ind(k)+1):(ind(k+1)-1));
			%
			if format_s(k)==1,
				eval(['s',int2str(k),'=champ;'],'1;');
			else
				indi=find(all((champ==' ')')');
				champn=champ;
				champn(indi,:)=char(ones(length(indi),1)*setstrl('nan',size(champn,2)));
				val=str2num(champn);
				if isempty(val),
					val=champ;
				end
				eval(['s',int2str(k),'=val;'],'1;');
			end
		end
	elseif ~isempty(findstr(txt(1,:),'Name:')),
		% 
		% recuperation des donnees sur help nom_table
		% separation du texte
		%
		texte1=txt(1:5,:);
		texte2=txt(8:size(txt,1),:);
		%
		% traitement du premier texte
		%
		nb1=size(texte1,1);
		for k=1:nb1,
			ind=min(find(texte1(k,:)==':'));
			tx1=texte1(k,1:(ind-1));
			tx2=texte1(k,(ind+1):size(texte1,2));
			eval(['s',int2str(k),'=tx2;'],'1;');
			nom_para(k,:)=setstrl(tx1,32);
		end
		%
		% deuxieme partie du texte
		%
		ind=all(texte2==' ');
		noml=texte2(1,:);
		noml_mem=noml;
		texte2=texte2(2:size(texte2,1),:);
		noml(noml==' ')=char(abs('_')*ones(1,sum(sum(noml==' '))));
		ftx=char(32*ind+noml.*(~ind));
		champ=[];
		tx=[' ',ftx,' '];
		tx(findstr(tx,'  '))=[];
		tx(findstr(tx,'  '))=[];
		ind=find(tx==' ');
		tx(tx=='_')=char(abs(' ')*ones(1,sum(sum(tx=='_'))));
		for k=1:(length(ind)-1)
			champ(k,:)=setstrl(tx((ind(k)+1):ind(k+1)),32);
		end
		nom_para=str2mat(nom_para,champ);
		nb_para=size(nom_para,1);
		lcm=size(champ,1);	
		%
		% recherche des indices
		%
		ind=[];
		for k=1:lcm,
			chh=champ(k,:);
			chh((1+max(find(chh~=' '))):length(chh))=[];
			ind(k)=findstr(noml_mem,chh);
		end
		ind(1)=1;
		ind(lcm+1)=size(texte2,2)+1;
		%
		% boucle sur les champs
		%
		for k=1:min(lcm,maxnbch),
			%
			% la variable au format texte
			%
			champ=texte2(:,ind(k):(ind(k+1)-1));
			indi=find(all((champ==' ')')');
			champn=champ;
			champn(indi,:)=char(ones(length(indi),1)*setstrl('nan',size(champn,2)));
			val=str2num(champn);
			if isempty(val),
				val=champ;
			end
			eval(['s',int2str(k+nb1),'=val;'],'1;');
		end		
	elseif ~isempty(findstr(txt(1,:),'Name')),
		%
		% cas du help sur une base
		%
		champ=[];
		tx=[' ',txt(1,:),' '];
		tx(findstr(tx,'  '))=[];
		tx(findstr(tx,'  '))=[];
		ind=find(tx==' ');
		for k=1:(length(ind)-1)
			champ(k,:)=setstrl(tx((ind(k)+1):(ind(k+1)-1)),32);
		end
		lcm=size(champ,1);
		nb_para=lcm;
		nom_para=champ;
		%
		% recherche des indices
		%
		ind=[];
		for k=1:lcm,
			chh=champ(k,:);
			chh(chh==' ')=[];
			ind(k)=findstr(txt(1,:),chh);
		end
		ind(1)=1;
		ind(lcm+1)=size(txt,2)+1;
		%
		% boucle sur les champs
		%
		for k=1:min(maxnbch,lcm),
			%
			% la variable au format texte
			%
			champ=txt(2:size(txt,1),ind(k):(ind(k+1)-1));
			indi=find(all((champ==' ')')');
			champn=champ;
			champn(indi,:)=char(ones(length(indi),1)*setstrl('nan',size(champn,2)));
			val=str2num(champn);
			if isempty(val),
				val=champ;
			end
			%
			% mise en forme
			%
			eval(['s',int2str(k),'=val;'],'1;');
		end	
	else
		disp('decodage impossible');
	end
