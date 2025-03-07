function cgclittcv(nom_data)

% initialisation
FF              =[];		
FP              =[];	 		
FPF             =[];			
LASTERROR       =[];	  
ans             =[];			
cr_exe          =[];		  
data            =[];	 	
dc              =[];			
ind             =[];			
liste           =[];		  
nom_exec        =[];	  
nomc            =[];	 	
noms            =[];	 	
numchoc         =[];	   
passwd          =[];		  
pid             =[];			
racine          =[];		  
user            =[];	  	
data            = {};

% chargement
[FP,FF] = fileparts(nom_data);
FPF     = fullfile(FP,[FF,'.log']);
diary(FPF);
load(nom_data)


% structure vide
dc.data = [];
dc.nom  = '';
dc.err  = NaN;

while ~isempty(liste)
	% nouveau signal
	nomc  =liste{end};
	liste(end) =[];
	% acces a la donnee
	fprintf('\n');
        [datac,err]   = mds_ftudata(numchoc,nomc);
        fprintf('\n');
        dc.data = datac;
        dc.nom  = nomc;
        dc.err  = err;
        % ecriture dans la structure de cellule
        data{end+1} = dc;
   
        % fin de la boucle
        save(nom_data)
        diary off
        diary on
end

