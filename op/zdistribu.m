% calcul distribue avec crnos
% ZDISTRIBU fonction feval distribuee
%------------------------------------------------------------------------------- 
% file :  zdistribu.m  ->   zdistribu
% 
% 
% function Matlab 7 : 
% 
% cette fonction permet d'appeler des fonction matlab en distribuant le calcul 
% sur differentes machines via le systeme de batch. L'appel est similaire a feval.
%
%
% syntaxe :  
%
%  appel d'une fonction distribue
%   
%    id = zdistribu('nom_fonction',entree1, ...,entreen);
%
%  attente de la fin des jobs :
%   
%   zdistribu;
%
%  lecture resultat :
%
%   [sortie1,...,sortien] = zdistribu(id);
%
%  arreter un job :
%
%    zdistribu(id,'KILL')
% 
%  arreter tous les jobs :
%    
%    zdistribu(0,'KILL')
%
%  information sur un job :
%
%   info =  zdistribu(id,'INFO')
%
%  information sur tous les jobs :
%
%   info =  zdistribu(0,'INFO')
%
%  fin de tous les jobs et reset du systeme :
%
%   zdistribu(0,'RESET')
% 
% input : 
%
%      'nom_fonction' :      nom de la fonction 
%      entree1..entreen :    variables d'entrees de la fonction
%      id  :                 numero d'identification du job
%  
% sorties :
%	
%      id  :                 numero d'identification du job
%      sortie1..sortien :    sorties de la fonction
%      info :                structure d'informations du job
%
% test :
%	id= zdistribu('testdistri',pi/3);  
%       zdistribu;
%       [e1,e2,e3] = zdistribu(id);
%
% fichiers a configurer sous architecture :
%
%  distribu_etat.sh  ->  retourne l'etat d'un job 
%  distribu_kill.sh  ->  arrete un job 
%  distribu_batch.sh -> lance un job
%  distribu.exe      -> execute un matlab distribue
%
% autres fonctions utiles :
%  testdistri.m      -> fonction de test
%  
% remarque : la session courante est sauve dans $HOME/.zdistribu.
% ce mecanisme ne permet uniquement a une seule session pour un meme user
% de redemarrer apres etre sortie de matlab.
%  
% function wrote by J-F Artaud, tel 62-15  
% version  3.3  du  13/04/2007
% CVS version  
%-------------------------------------------------------------------------------  
%  
function varargout = zdistribu(entree1,varargin)

% recuperation variable environement
if isdir(fullfile(getenv('HOME'),'zineb','data'));
    workdir = fullfile(getenv('HOME'),'zineb','data');
else
    workdir  = pwd;
end

% donnees
if isappdata(0,'BATCH_DISTRIBUE_CRONOS')
	ba = getappdata(0','BATCH_DISTRIBUE_CRONOS');
elseif exist(fullfile(getenv('HOME'),'.zdistribu.mat'),'file');
	ba =struct([]);
	warning off
		load(fullfile(getenv('HOME'),'.zdistribu.mat'));
	warning on
else
	ba =struct([]);
end
% commande de lancement et de controle des jobs
cmd.run   = 'architecture/distribu_batch.sh';
cmd.state = 'architecture/distribu_etat.sh';
cmd.kill  = 'architecture/distribu_kill.sh';

% root de cronos
root        = getappdata(0,'root');
cmd.run     = fullfile(root,cmd.run);
cmd.state   = fullfile(root,cmd.state);
cmd.kill    = fullfile(root,cmd.kill);

% pas de sortie par defaut
varargout = {};

% le mode depend des entrees
if nargin == 0
	% test des etats et attente
	attente = 1;
	while	(attente == 1)
		attente = 0;
		% mise a jour des etat
		for k = 1:length(ba)
			switch ba(k).etat 
			case 'running'
				[s,t] = unix(sprintf('%s %d',cmd.state,ba(k).jobid));
				if s ~= 0
					ba(k).etat ='done';	
					fprintf('ZDISTRIBU: job %d terminated (%d)\n',k,ba(k).jobid);
				else
					attente = 1;
					pause(0.1);
				end
			end
		end
	end

elseif isstr(entree1)
	% id du job
	current = length(ba) +1;	
	% rempplissage de la structure
	% commande
	ba(current).commande = entree1;
	% donnees de la commande
	%ba(current).entrees   = varargin;
	% chemin matlab
	%ba(current).path     = path;
	% parametre de l'application
	%ba(current).appdata  = getappdata(0);
	% repertoire d'execution
	ba(current).dir      = pwd;
	% fichier journal
	ba(current).diaryfile   = fullfile(workdir,sprintf('diary_distribu_%d.txt',current));
	% fichier des resultats
	ba(current).outputdata   = fullfile(workdir,sprintf('output_distribu_%d.mat',current));
	% fichier de lencement du batch
	ba(current).batchfile   = fullfile(workdir,sprintf('batch_distribu_%d',current));
	% fichier de lencement du batch
	ba(current).batchout   = fullfile(workdir,sprintf('zd%d.o*',current));
	% etat du run
	ba(current).etat  = 'prepare';
	% sortie en erreur
	ba(current).etat_sortie  = 1;
	% message d'erreur
	ba(current).error  = '';
	% texte ecrit pendant le run
	%ba(current).diary  = '';
	% resultat
	%ba(current).resultat  = {};
	
	% appel equivalent a feval
	% sauvegarde des donnees
	ba(current).inputdata   = fullfile(workdir,sprintf('input_distribu_%d.mat',current));
	data                    = ba(current);
	data.entrees            = varargin;
	data.path               = path;
	data.appdata            = getappdata(0);
	data.resultat           = {};
	save(ba(current).inputdata,'data');
	% lancement du job
	script = fullfile(root,'architecture','distribu.exe');
	txt =   sprintf('%s %s zd%d %s %s %s %s', ...
	        cmd.run,workdir,current,getenv('USER'),script,ba(current).inputdata,ba(current).batchfile);
	[s,t] = unix(txt);
	if s~= 0
	  	error(sprintf('ZDISTRIBU : error launching job %d in batch :\n %s',entree1,t));				
	else
		[v,r] = strtok(t);
		[v,r] = strtok(r);
		[v,r] = strtok(r);
		ba(current).jobid = str2num(v);	
		ba(current).etat  = 'running';

	end
	varargout{1} = current;
elseif nargin == 1
	% lecture du resultat pour la reference donnees
	if length(ba) >= entree1
		switch ba(entree1).etat 
		case 'running'
			[s,t] = unix(sprintf('%s %d',cmd.state,ba(entree1).jobid));
			if s ~= 0
				ba(entree1).etat ='done';	
			end
		end
		% test de fin
		switch ba(entree1).etat 
		case  'done'
			% lecture du resultat
			try 
				data = load(ba(entree1).outputdata);
			catch
				if exist(ba(entree1).diaryfile,'file')
					type(ba(entree1).diaryfile)
				end
				error(sprintf('ZDISTRIBU : error loading result job %d',entree1));			
			end
			% lecture du texte
			%ba(entree1).resultat = data.resultat;
			ba(entree1).etat_sortie = data.etat_sortie;
			ba(entree1).error = data.etat_error;
			
			if ba(entree1).etat_sortie == 1
				error(sprintf('ZDISTRIBU : error in "%s" :\n %s',  ...
				      ba(entree1).commande,ba(entree1).error));
			end
			varargout = data.resultat;
			% destruction du resultat
			warning off
			delete(ba(entree1).diaryfile);
			delete(ba(entree1).inputdata);
			delete(ba(entree1).outputdata);
			delete(ba(entree1).batchfile);
			delete(ba(entree1).batchout);
			warning on
			ba(entree1).etat = 'delete';
		case 'delete'
			warning('ZDISTRIBU : job %d results had already been read',entree1);
			[varargout{1:nargout}] = deal([]);
		otherwise
			warning('ZDISTRIBU : job %d is not finish',entree1);
			[varargout{1:nargout}] = deal([]);
		end
	end
else
	if entree1 == 0
		% test des etats et attente
		attente = 1;
		while	(attente == 1)
			attente = 0;
			% mise a jour des etat
			for k = 1:length(ba)
				switch ba(k).etat 
				case 'running'
					[s,t] = unix(sprintf('%s %d',cmd.state,ba(k).jobid));
					if s ~= 0
						ba(k).etat ='done';	
					end
				end
			end
		end
	else
		switch ba(entree1).etat 
		case 'running'
			[s,t] = unix(sprintf('%s %d',cmd.state,ba(entree1).jobid));
			if s ~= 0
				ba(entree1).etat ='done';	
			end
		end
	end
	
	switch varargin{1}
	case 'KILL'
		if entree1 == 0
			for k = 1:length(ba)
				switch ba(k).etat 
				case 'running';
					[s,t] = unix(sprintf('% %d',cmd.kill,ba(k).jobid));
					% destruction du resultat
					warning off
					delete(ba(k).diaryfile);
					delete(ba(k).inputdata);
					delete(ba(k).outputdata);
					delete(ba(k).batchfile);
					delete(ba(k).batchout);
					warning on
					ba(k).jobid = NaN;
					ba(k).etat ='delete';
				case 'delete'
					%rien
				otherwise
					warning off
					delete(ba(k).diaryfile);
					delete(ba(k).inputdata);
					delete(ba(k).outputdata);
					delete(ba(k).batchfile);
					delete(ba(k).batchout);
					warning on
					ba(k).jobid = NaN;
					ba(k).etat ='delete';					
				end	
			end	
		elseif length(ba) >= entree1 
			% fin du job					
			switch ba(entree1).etat 
			case 'running'
				[s,t] = unix(sprintf('% %d',cmd.kill,ba(entree1).jobid));
				if s == 0
					ba(entree1).jobid = NaN;
					% destruction des fichiers 
					warning off
					delete(ba(entree1).diaryfile);
					delete(ba(entree1).inputdata);
					delete(ba(entree1).outputdata);
					delete(ba(entree1).batchfile);	 
					delete(ba(entree1).batchout);
					warning on
					ba(entree1).etat ='delete';							
				else
					error(sprintf('ZDISTRIBU : unable to delete job %',ba(entree1).jobid));
				end
			case 'delete'
				%rien
			otherwise
				ba(entree1).jobid = NaN;
				% destruction des fichiers 
				warning off
				delete(ba(entree1).diaryfile);
				delete(ba(entree1).inputdata);
				delete(ba(entree1).outputdata);
				delete(ba(entree1).batchfile);	 
				delete(ba(entree1).batchout);
				warning on
				ba(entree1).etat ='delete';
										
			end
		end
	case 'INFO'
		if entree1 == 0
			varargout{1} =	ba;
		else
			varargout{1} =	ba(entree1);		
		end	
	case 'RESET'
		zdistribu(0,'KILL')
		ba = struct([]);
	end
	
end

setappdata(0,'BATCH_DISTRIBUE_CRONOS',ba);	
save(fullfile(getenv('HOME'),'.zdistribu.mat'),'ba');
