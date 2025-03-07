% script de chargement des donnees
if ~exist('data','var') | ~exist('param','var')
     zuiload;
end
if exist('jeux1','var')
        if isfield(jeux1,'param')
	  if fix(param.from.shot.num) ~= fix(jeux1.param.from.shot.num)
		clear jeux1;
	  end
	else
           clear jeux1
        end
	
end

if ~exist('jeux1','var') 
      mask = strcat('*',int2str(fix(param.from.shot.num)),'*.mat*');
      [file,pathf]=uigetfile(mask,'reference data ?');
      drawnow
      if ischar(file)
            [lfile,rmfile,cr]=zgzip(strcat(pathf,file),'uncompress');
            jeux1 = load(lfile);
            jeux1.data=zreduit(jeux1.param,jeux1.data,'uncompact');
            if ~isempty(rmfile)
            	[voids,voidt]=unix(['rm -f ',rmfile,' >& /dev/null']);
            end
            %jeux1 = load(strcat(pathf,file));
            t1 = jeux1.data.gene.temps;  
     else
            jeux1=[];
            t1 =[];  
      end
end

switch param.from.machine

case 'TS'
	if exist('bile','var') | exist('sbile','var')
	        if ~isempty(bile.ip)
		  if str2num(bile.info.bile.choc) ~= fix(param.from.shot.num)
			clear bile sbile
		  end
		end
	end
	
	if ~exist('bile','var') | ~exist('sbile','var')
  	   [bile,sbile] =cgcgettrait(param.from.shot.num,'tprof');
	end
	% 1- preparation du changement de rho
	vtin = ones(size(bile.times));
	vein = ones(size(bile.rhofit));
	rhogin =  vtin * bile.rhofit;
	shiftin = (bile.d0 * vein) .* ( 1 - rhogin .^ (bile.piqd * vein));
	delta = ((bile.rmaj * vein) + shiftin) .^ 2 - (bile.amin * vein).^2 .* rhogin.^2;
	phiin = 2 .* pi .* (bile.rmaj * vein) .* (bile.btor * vein) .* ...
       				 ((bile.rmaj * vein) + shiftin - sqrt(delta));
        
	% coordonnee de flux dans la base temps d'origine        
	rhoin = sqrt(phiin ./ pi ./ (bile.btor * vein));   % coordonnee sqrt(Phi/pi/B0) pour les donnees bile
	rhomaxin = max(rhoin')';   % rhomax dans bile
	xbile = rhoin ./ (rhomaxin * ones(1,size(rhoin,2)));
        times = bile.times;

        % flux de neutron
        if ~exist('ndd_exp','var')
           [ndd,tndd,void1,void2] = tsbase(fix(param.from.shot.num),'gfluntn');
           if ~isempty(ndd)
              if param.compo.z(1) == 1
                   k = 1;
              else
                   k = 2;
              end               
              ndd_exp  = 10 .^ ndd(:,k);
              tndd_exp = tndd(:,k);
           end
        end

        % tiprof
        if ~exist('tiprof','var')
           numchoc = param.from.shot.num;
           file  = sprintf('/usr/drfc/cgc/matlab5/tcron/TS/%d/Ti%d.mat',fix(numchoc),fix(numchoc));

           if exist(file,'file')
              tiprof=load(file);
           else
              disp('fichier de donnees tiprof indisponible')
              tiprof=[];
           end
        
       end
case 'JET'
        if exist('jettemp','var')
		  if isfield(jettemp,'numchoc')
	       if param.from.shot.num ~= jettemp.numchoc
	         clear jettemp jetprof jettethom jetticx jetsh jetteshtom jettransp jetefit
	       end
		  end
	end
	if strcmp(getenv('HOME'),'/usr/drfc/cgc')
	  chemin =  strcat(getenv('HOME'),'/matlab5/tcron/JET/data/');
	else
	  chemin =  strcat(getenv('HOME'),'/zineb/data/JET/');
	  if ~exist(chemin,'dir')
	    disp(['the directory ',chemin,'does not exist!'])
		 chemin = input('new  directory ','s')
	  end
	end
	if ~exist('jettemp','var')
	   nomfichier = sprintf('%s%d/temp%d',chemin, ...
	                        param.from.shot.num,param.from.shot.num);
           jettemp = load(nomfichier);
        end
	if ~exist('jetprof','var')
	   nomfichier = sprintf('%s%d/prof%d',chemin,...
	                        param.from.shot.num,param.from.shot.num);
           jetprof = load(nomfichier);
        end
	if ~exist('jettethom','var')
	   nomfichier = sprintf('%s%d/tex%d',chemin, ...
	                        param.from.shot.num,param.from.shot.num);
	   if exist(nomfichier,'file')			
             jettethom = load(nomfichier);
	   else
	     jettethom = [];
	   end
	end
	if ~exist('jetticx','var')
	   nomfichier = sprintf('%s%d/tix%d',chemin, ...
	                        param.from.shot.num,param.from.shot.num);
	   if exist(nomfichier,'file')			
             jetticx = load(nomfichier);
	   else
	     jetticx = [];
	   end
	end
	if ~exist('jetsh','var')
	   nomfichier = sprintf('%s%d/teshx%d',chemin, ...
	                        param.from.shot.num,param.from.shot.num);
	   if exist(nomfichier,'file')			 
             jetsh = load(nomfichier);
	   else
	     jetsh = [];
	   end
	end
	
	if ~exist('jetteshthom','var')
	   nomfichier = sprintf('%s%d/teshthx%d',chemin, ...
	                        param.from.shot.num,param.from.shot.num);
	   if exist(nomfichier,'file')			 				
             jetteshthom = load(nomfichier);
	   else
	     jetteshthom = [];
	   end
	end		
	if ~exist('jettransp','var')
	       nomfichier = sprintf('%s%d/transx%d',chemin, ...
	                              param.from.shot.num,param.from.shot.num);
	       if exist(strcat(nomfichier,'.mat'),'file')
	          jettransp = load(nomfichier);
	       else
	          jettransp = [];
	       end
	end
	if ~exist('jetefit','var')
	       nomfichier = sprintf('%s%d/efit%d',chemin, ...
	                              param.from.shot.num,param.from.shot.num);
	       if exist(strcat(nomfichier,'.mat'),'file')
	          jetefit = load(nomfichier);
	       else
	          jetefit = [];
	       end
	end
case 'DIIID'
     home                           = getenv('HOME');
     racine                         = cat(2,home,'/zineb/data/diiid/');
     directory                      = [racine,int2str(param.from.shot.num)];
     if ~exist(directory)
     
       directory = input('directory name for the DIIID input files :','s');
       directory = fullfile(directory,int2str(param.from.shot.num));
     
     end
     filetemp = fullfile(directory,'diiidtemp.mat');
     fileprof = [directory,'/diiidprof.mat'];
     filediag = [directory,'/diiiddiag.mat'];
     eval(['load ',filetemp])
     eval(['load ',fileprof])
case 'HL2A'
     home                           = getenv('HOME');
     racine                         = cat(2,home,'/zineb/data/hl2a/');
     directory                      = [racine,int2str(param.from.shot.num)];
     filetemp = [directory,'/hl2atemp.mat'];
     fileprof = [directory,'/hl2aprof.mat'];
     fileeq = [directory,'/hl2aeq.mat'];
     load(filetemp)
     load(fileprof)
     load(fileeq)

otherwise
        disp('No input data for this tokamak')
end

x =param.gene.x;
t =data.gene.temps;
if ~exist('t1','var')
  t1 = [];
end
% compatibilite ascendante
if ~isempty(t1)
   if ~isfield(jeux1.data.gene,'ihyb')
      jeux1.data.gene.ihyb = NaN .* jeux1.data.gene.temps;
      disp('fichier ancien')
   end
end

if ~isempty(t1)
   if ~isfield(jeux1.data.gene,'iidn')
      jeux1.data.gene.iidn = NaN .* jeux1.data.gene.temps;
      disp('fichier ancien')
   end
end
if ~isfield(data.gene,'ihyb')
   data.gene.ihyb = NaN .* data.gene.temps;
   disp('old file')
end

if ~isfield(data.gene,'iidn')
   data.gene.iidn = NaN .* data.gene.temps;
   disp('old file')
end

