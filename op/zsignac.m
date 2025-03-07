% ZSIGNAC change le nombre de points radiaux d'un jeu de donnees CRONOS
%-----------------------------------------------------------------------
%
% fichier zsignac.m
%
% fonction Matlab 5 :
%	Cette fonction change le nombre de point radiaux d'un jeu de donnees CRONOS.
% 
% syntaxe
%        [data_new,param_new] = zsignac(data,param,option);
% entrees :
%	param    = structure de parametress de cronos
%       data     = structure de donnees de cronos
%       option   = structure des option de la fonction (cf partie autodeclarante)
%
% sorties :
%       data_new   = nouvelle structure de donnees de cronos (avec le nombre de points radiaux souhaites)
%       param_new  = nouvelle structure de parametress de cronos
%
% fonction ecrite par J-F Artaud, poste 62-15
% version  2.1 , du 17/06/2003.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function [data_new,param_new] = zsignac(data,param,option)

if nargin == 0
	langue                  = getappdata(0,'langue_cronos');
	
	valeur.nbrho   = 101;                        % 
	type.nbrho     = 'integer';     % type de la donnnee
	borne.nbrho    = [21,1001];         % bornes
	defaut.nbrho   = 101;           % valeur par defaut
        info.nbrho   = 'nombre de points radiaux souhaites [21,1001], impair';

	valeur.methode    = 'auto';                        % 
	type.methode     = 'string';     % type de la donnnee
	borne.methode    =  {'auto','nearest','linear','spline'};       % bornes
	defaut.methode   = 'auto';           % valeur par defaut
        info.methode   = 'methode d''interpolation (cf. interp1)';
	sortie.description = 'Change le nombre de points radiaux d''un jeu de donnees cronos';   % description (une ligne) de la fonction


        if strcmp(langue,'anglais')	
          info.nbrho   = 'number of radial points [21,1001], even';
          info.methode   = 'Interpolation method (i.e. interp1)';
	sortie.description = 'modify the radial coordinate x of the profiles';   % description (une ligne) de la fonction
        end
	
	
	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	
	sortie.help = 'help zsignac';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	data_new = sortie;
	
	return

end


% compatibilite
nbrho_new = option.nbrho;
if isempty(nbrho_new)
   nbrho_new = 21;
end
nbrho_new = 2.* fix(nbrho_new/2) +1;

methode   = option.methode;
if strcmp(methode,'auto')
   methode ='';
end
if isempty(methode)
   if nbrho_new > length(param.gene.x)
      methode = 'spline';
   elseif nbrho_new < fix(length(param.gene.x)/3);
      methode = 'nearest';
   else
      methode = 'linear';
   end
end

% nouvelles donnees
[cr,data_new,param_new]=zinit('',data.gene.temps,nbrho_new,'',[],[], ...
                                 param.nombre.fci,param.nombre.fce,param.nombre.idn, ...
                                 param.nombre.hyb,param.nombre.glacon,param.gene.nbg);
if cr ~=0
	return
end

% connexion des modules externes
[cr,data_new,param_new] = zconnexion(data_new,param_new);
if cr ~=0
	return
end

% nouvelle structure param   
param_mem  = param_new;                              
param_new  = param;
param_new.gene.nbrho =  nbrho_new; 
param_new.gene.x =  param_mem.gene.x; 
param_new.asser  =  param_mem.asser; 
param_new.memoire  =  param_mem.memoire; 
param_new.batch  =  param_mem.batch; 



                               
% boucle sur les champs de data
disp('Data');
% liste des champs la structures
champ = fieldnames(data);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('data.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while(~isempty(champ))
    %premier champ de la liste
    champc=champ{1};
    champ(1)=[];
    eval(strcat('test=isstruct(',champc,');'));
    if test
    	% cas d'une sous structure -> ajout de champs
    	eval(strcat('champnew=fieldnames(',champc,');'));   	
    	for k=1:length(champnew)
    		% nom complet pour acceder a la variable
	        champnew{k}=strcat(champc,'.',champnew{k});
	   end
	   % ajout a la liste des champs
	   if isempty(champ)
	      champ =champnew;
	   else
	      champ=cat(1,champ,champnew);
	   end
    else
      % juste pour les tests
      disp(champc);
      % valeur
      val = eval(champc);
      % nouveau nom
      new = strrep(champc,'data.','data_new.');
      % interpolation
      if ~isempty(findstr(champc,'rhoRZ'))
           % exception, on ne fait rien
      elseif length(size(val))> 2
         if ~isempty(findstr(champc,'.equi.'))
            % donnees grille equilibre on ne fait rien
         else
            valm = val;
            val  = NaN .* ones(length(data_new.gene.temps),length(param_new.gene.x), size(valm,3));
            for k = 1:size(val,3)
                  vv = squeeze(valm(:,:,k));
                  val(:,:,k) = interp1(param.gene.x',vv',param_new.gene.x',methode)';
            end
         end
      elseif size(val,2) > 1
         if length(param.gene.x) == size(val,2)
            val = interp1(param.gene.x',val',param_new.gene.x',methode)';
         end
      else
         % rien;
      end
      
      eval(sprintf('%s = val;',new));
      
    end
end
