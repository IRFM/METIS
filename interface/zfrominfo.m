% ZFROMINFO affiche liste info de param.from.shot.info et param.from.creation.info
%----------------------------------------------------------------
% fichier zfrominfo.m 
%
% fonction Matlab 5 :
%	Cette fonction permet d'afficher les donnees "param.from.shot.info
%	et "param.from.creation.info" de zineb.
%	Les donnees doivent etre dans l'espace de travail de base.
%	Elles sont contenue dans la structure.param
% 
% syntaxe  :
%	tmp = zfrominfo(ind,{param});
%
% entrees
%	ind : =1   on traite  param.from.shot.info 
%	      =2   on traite param.from.creation.info
%
% sorties
%	tmp : nom du fichier temporaire cr�
%
% fonction �rite par C. Passeron , poste 61-19
% version  1.7  du  29/09/2001  
%  
% liste des modifications : 
%	* 12/07/2001  -> le fichier temporaire sous le repertoire /tmp est cree 
%                   en fonction du pid du process Matlab     
%
%  * 18/07/2002 -> ajout du 2ieme argument
%--------------------------------------------------------------
%
function tmp = zfrominfo(ind,param)


if nargin<1
	disp('syntaxe:  tmp = zfrominfo(ind,{param})')
	return
end
if ind~=1 & ind~=2
	disp('first input must be 1 or 2') ;
	return
end
if nargin < 2
   param =[];
end

% Ouverture d'un fichier temporaire
%tmp = tempname ;
[pid,pid_sess,user]=getidprocess ;
tmp = strcat('/tmp/info',num2str(pid));
[fid,mess] = fopen(tmp,'w') ;
if ~isempty(mess)
	herror = errordlg('error opening file','Warning') ;
	zwaitfor(herror) ;
	zuicloseone(hform) ;
end

switch ind 
case {1}

	% liste des donnees
	if isempty(param)
		try 
			param.from.shot.info = evalin('base','param.from.shot.info') ;
		catch
			param.from.shot.info=[] ;
		end
	end
	if isempty(param.from.shot.info)  
		herror = errordlg('no data in "param.from" ','Warning') ;
		zwaitfor(herror) ;
		zuicloseone(hform) ;
		return
	end

	if ~isstruct(param.from.shot.info)
		return
	end
	 
	fprintf(fid,' param.from.shot.info \n -------------------- \n') ;

	% liste des champs la structures
	champ = fieldnames(param.from.shot.info) ;
	nom   = fieldnames(param.from.shot.info) ; 

	for k=1:length(nom)
		% nom complet pour acceder a la variable
		champ{k} = strcat('param.from.shot.info.',champ{k}) ;
	end

	test = 0 ;
	while (~isempty(champ))
		fprintf(fid,'%s \n',nom{1}) ;
		nom(1)=[] ;
	
		%premier champ de la liste
		champc=champ{1} ;
		champ(1)=[] ;
		eval(strcat('test=isstruct(',champc,');')) ;
		if test
	   		% cas d'une sous structure -> ajout de champs
	   		eval(strcat('champnew=fieldnames(',champc,');')) ;
		
	   		for k=1:length(champnew)
				nom2 = champnew{k} ;
	   			% nom complet pour acceder a la variable
	   			champnew{k}=strcat(champc,'.',champnew{k}) ;
			
				eval(strcat('test=isstruct(',champnew{k},');')) ;
				if test
	   				% cas d'une sous structure -> ajout de champs
	   				eval(strcat('champnew2=fieldnames(',champnew{k},') ;')) ;
				
			   		for l=1:length(champnew2)
						nom3 = champnew2{l} ;
						% nom complet pour acceder a la variable
			   			champnew2{l}=strcat(champnew{k},'.',champnew2{l}) ;
					
						eval(strcat('test=isstruct(',champnew2{l},');')) ;
						if test
	   						% cas d'une sous structure -> ajout de champs
	   						eval(strcat('champnew3=fieldnames(',champnew2{l},');')) ;
						
			   				for m=1:length(champnew3)
								nom4 = champnew3{m} ;
								% nom complet pour acceder a la variable
			   					champnew3{m}=strcat(champnew2{l},'.',champnew3{m}) ;
							
								fprintf(fid,'\t %s.%s.%s \t: %s \n',nom2,nom3,nom4,eval(champnew3{m})) ;
							end
						else
					
							fprintf(fid,'\t %s.%s \t: %s \n',nom2,nom3,eval(champnew2{l})) ;
						end
					end
				else
					fprintf(fid,'\t %s \t: %s \n',nom2,eval(champnew{k})) ;
				end
				
	   		end
	   
		end
	end         

case {2}
	%disp('processing param.from.creation.info')

	% liste des donnees
	if isempty(param)
		try 
			param.from.creation.info = evalin('base','param.from.creation.info') ;
		catch
			param.from.creation.info=[] ;
		end
	end
	if isempty(param.from.creation.info)  
		herror = errordlg('no data in "param.from" ','Warning') ;
		zwaitfor(herror) ;
		zuicloseone(hform) ;
		return
	end

	if ~isstruct(param.from.creation.info)
		return
	end
	 
	fprintf(fid,' param.from.creation.info \n -------------------- \n') ;

	% liste des champs la structures
	champ = fieldnames(param.from.creation.info) ;
	nom   = fieldnames(param.from.creation.info) ; 

	for k=1:length(nom)
		% nom complet pour acceder a la variable
		champ{k} = strcat('param.from.creation.info.',champ{k}) ;
	end

	test = 0 ;
	while (~isempty(champ))
		fprintf(fid,'%s \n',nom{1}) ;
		nom(1)=[] ;
	
		%premier champ de la liste
		champc=champ{1} ;
		champ(1)=[] ;
		eval(strcat('test=isstruct(',champc,');')) ;
		if test
	   		% cas d'une sous structure -> ajout de champs
	   		eval(strcat('champnew=fieldnames(',champc,');')) ;
		
	   		for k=1:length(champnew)
				nom2 = champnew{k} ;
	   			% nom complet pour acceder a la variable
	   			champnew{k}=strcat(champc,'.',champnew{k}) ;
			
				eval(strcat('test=isstruct(',champnew{k},');')) ;
				if test
	   				% cas d'une sous structure -> ajout de champs
	   				eval(strcat('champnew2=fieldnames(',champnew{k},') ;')) ;
				
			   		for l=1:length(champnew2)
						nom3 = champnew2{l} ;
						% nom complet pour acceder a la variable
			   			champnew2{l}=strcat(champnew{k},'.',champnew2{l}) ;
					
						eval(strcat('test=isstruct(',champnew2{l},');')) ;
						if test
	   						% cas d'une sous structure -> ajout de champs
	   						eval(strcat('champnew3=fieldnames(',champnew2{l},');')) ;
						
			   				for m=1:length(champnew3)
								nom4 = champnew3{m} ;
								% nom complet pour acceder a la variable
			   					champnew3{m}=strcat(champnew2{l},'.',champnew3{m}) ;
							
								fprintf(fid,'\t %s.%s.%s \t: %s \n',nom2,nom3,nom4,eval(champnew3{m})) ;
							end
						else
					
							fprintf(fid,'\t %s.%s \t: %s \n',nom2,nom3,eval(champnew2{l})) ;
						end
					end
				else
					fprintf(fid,'\t %s \t: %s \n',nom2,eval(champnew{k})) ;
				end
				
	   		end
	   
		end
	end         

otherwise
	disp('first input must be 1 or 2') ;
	return
end
fclose(fid) ;

