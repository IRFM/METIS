% ZINEB1TEST effectue le test de un module de CRONOS ou de certaines fonction de CRSONOS
%----------------------------------------------------------------------------------------
% fichier zineb1test.m ->  zineb1test,locmergeparam,zsetfield
%
%
% fonction Matlab 7 :
%
% Cette fonction permet d'appeler les modules de CRONOS ou de certaine focntion de CRONOS et
% de comparer le resultat avec celui obtenu lors d'un run precedent.
% La fonction affiche les differences entre les differentes donnees (format texte et/ou graphique)
%
% syntaxe  :
%
%   pour effectuer les tests :
%
%      [datak,param] = zineb1test(data,param,temps,module,{liste de pair d'options});
%
%   pour avoir la liste des tests possibles :
%
%  liste = zineb1test
%
%
% entrees :
%
%     data   =  structure data de CRONOS (au moins un temps pour les tests de modules, sinon au moins 2 temps)
%
%     param  = structure param de CRONOS.
%
%     temps  = temps chosi pour le test (s)
%
%     module = mot clef pour le test d'un module (liste = zineb1test)
%
% options :
%
%     ,'reprise',nombre  = force le mode reprise de CRONOS
%     ,'option',structure option du module = change le parametrage du module
%     ,'nomf',nom du module = change le module connecte a CRONOS
%     ,'plotonoff, 1  = cree les fichiers .png des plot de comparaison des donnees
%     ,'plotonoff, 0  = sorties textes seulement
%     ,'updateparam',0 = pas de mise a jour automatique des parametres d'un module
%     ,'updateparam',1 = pas mise a jour automatique des parametres d'un module
%     ,'debug',0       = pas d'arret si plantage d'un module
%     ,'debug',1       = arret lors du plantage d'un module
%     ,'timestep',delta_t = pas de temps optionnel pour l'option onestep
%     ,'clean',0       = n'efface pas les fichiers cree par les modules
%     ,'clean',1       = efface les fichiers cree par les modules
%     ,'tolerance',1e-3  = tolerance sur le resultat
%
%
%
% sorties :
%
%     datak    =   structure data a 1 temps contenant les donnees modifiees .
%
%     param = structure param modifiee.
%
% fonction ecrite par J-F Artaud , poste 62-15
% version 3.1, du 23/10/06.
%
%
% liste des modifications :  version CVS
%
%--------------------------------------------------------------
%
function [datak,param,tolerance] = zineb1test(data,param,temps,module,varargin)

% securite
clear global

% delaration de la liste des modules testables
if nargin == 0
        datak.equi      = 'test of equilibrium module';
        datak.neo       = 'test of neoclassical module';
        datak.impur     = 'test of impurities module';
        datak.rip       = 'test of ripple module';
        datak.mhd_dds      = 'test of dds reconnexion module';
        datak.mhd_limite   = 'test of MHD reconnexion limite module';
        datak.mhd_elm      = 'test of ELMs module';
        datak.mhd_stab     = 'test of MHD stability analyze module';
        datak.fci    = 'test of ICRH module';
        datak.fce    = 'test of ECRH module';
        datak.hyb    = 'test of LHCD module';
        datak.idn    = 'test of NBI module';
        datak.n0                = 'test of cold edge neutrals and ionisation source module';
        datak.bord              = 'test of SOL module';
        datak.glacon            = 'test of pellet module';
        datak.fus               = 'test of fusion sources module';
        datak.cyclo     = 'test of synchrotron radiation module';
        datak.coefa     = 'test of transport coefficient module A';
        datak.coefb     = 'test of transport coefficient module B';
        datak.coefc     = 'test of transport coefficient module C';
        datak.coefd     = 'test of transport coefficient module D';
        datak.coefe     = 'test of transport coefficient module E';
        datak.coeff     = 'test of transport coefficient module F';
        %datak.plot             = 'not yet implemented';
        datak.asser             = 'test of feedback control module';
        %datak.post             = 'not yet implemented';
        datak.coef_all          = 'test of all connected transport coefficient modules';
        datak.sources           = 'test of all sources modules';
        datak.first             = 'test of initialisation phase of a CRONOS run';
        datak.onestep           = 'test of cronos on one time step';
        param                   = datak;
   return

end

% initialisation des sorties
datak =[];

%force l'interactif
setappdata(0,'inter',1);

% selon les entrees
% test du type
listoftype = fieldnames(zineb1test);
if isempty(strmatch(module,listoftype,'exact'))
   fprintf('type of module unknown (%s) !\n',module);
   return
end

% changement de syntaxe
module_st = strrep(module,'_','.');


% nomf,option,reprise
if nargin < 5
   varargin =[];
end
nomf    = '';
option  = [];
if isfield(param.gene,'reprise')
   reprise = param.gene.reprise;
else
   reprise = 0;
end
plotonoff = 0;
debug_onoff = 0;
updateparam = 0;
tolerance = 1e-3;
timestep = [];
clean    = 1;
dbclear all
if ~isempty(varargin)
   ind = strmatch('option',varargin(1:2:end-1),'exact');
   if ~isempty(ind)
      option = varargin{2 .* ind};
   end
   ind = strmatch('nomf',varargin(1:2:end-1),'exact');
   if ~isempty(ind)
      nomf = varargin{2 .* ind};
   end
   ind = strmatch('reprise',varargin(1:2:end-1),'exact');
   if ~isempty(ind)
      try
         reprise = eval(varargin{2 .* ind});
      catch
         reprise = varargin{2 .* ind};
      end
   end
   ind = strmatch('plotonoff',varargin(1:2:end-1),'exact');
   if ~isempty(ind)
      try
         plotonoff = eval(varargin{2 .* ind});
      catch
         plotonoff = varargin{2 .* ind};
      end
   end
   ind = strmatch('debug',varargin(1:2:end-1),'exact');
   if ~isempty(ind)
      try
         debug_onoff = eval(varargin{2 .* ind});
      catch
         debug_onoff = varargin{2 .* ind};
      end
   end
   ind = strmatch('updateparam',varargin(1:2:end-1),'exact');
   if ~isempty(ind)
      try
         updateparam = eval(varargin{2 .* ind});
      catch
         updateparam = varargin{2 .* ind};
      end
   end
   ind = strmatch('timestep',varargin(1:2:end-1),'exact');
   if ~isempty(ind)
      try
         timestep = eval(varargin{2 .* ind});
      catch
         timestep = varargin{2 .* ind};
      end
   end
   ind = strmatch('clean',varargin(1:2:end-1),'exact');
   if ~isempty(ind)
      try
         clean = eval(varargin{2 .* ind});
      catch
         clean = varargin{2 .* ind};
      end
   end
   ind = strmatch('tolerance',varargin(1:2:end-1),'exact');
   if ~isempty(ind)
      try
         tolerance = eval(varargin{2 .* ind});
      catch
         tolerance = varargin{2 .* ind};
      end
   end
end
% update
param.gene.reprise = reprise;


% poure la comparaison
memin    = param.memoire;
if isfield(param,'assermem');
   assermem = param.assermem;
else
   assermem = param.asser;
end
%
% positionnement correct des variables memoires en reprise
%
if (reprise == 1) | (reprise == 0)
   param.memoire.fci.data   = [];
   param.memoire.idn.data   = [];
   param.memoire.hyb.data   = [];
   param.memoire.fce.data   = [];
   param.memoire.impur.data = [];
   param.memoire.bord.data  = [];
   param.memoire.fus.data   = [];
   param.memoire.n0.data    = [];
   param.memoire.stab.data  = [];
   param.memoire.equi.data  = [];
   param.memoire.neo.data   = [];
   param.memoire.cyclo.data = [];
   param.asser              = [];
   param.asser.data         = [];
   param.asser.delais       = [];
   premiertemps = 1;
else
   premiertemps = 0;
end

% effacement du contexte global
setappdata(0,'equilibre_helena_init',[]);
setappdata(0,'equilibre_helena_b',[]);
setappdata(0,'equilibre_helena_df2',[]);
setappdata(0,'equilibre_helena_kkbig',[]);
setappdata(0,'equilibre_helena_ptotold',[]);
setappdata(0,'equilibre_helena_corj',[]);
setappdata(0,'sinbad_data',[]);
rmappdata(0,'equilibre_helena_init');
rmappdata(0,'equilibre_helena_b');
rmappdata(0,'equilibre_helena_df2');
rmappdata(0,'equilibre_helena_kkbig');
rmappdata(0,'equilibre_helena_ptotold');
rmappdata(0,'equilibre_helena_corj');
rmappdata(0,'sinbad_data');

% on force le mode verbose
param.gene.verbose = 1;
setappdata(0,'zverbose',1);

% separateur
disp('--------------------------------------------------')


% connexion du module si necessaire
if ~isempty(nomf)
   try
      nb = zgetfield(param.nombre,module);
   catch
      nb = 1;
   end
   param.fonction = zsetfield(param.fonction,module_st,nomf);
   %[cr,data,param] = zconnexion(data,param);
        [param,nbok,cr] =  locmergeparam(param,module_st,zgetfield(param.fonction,module_st),nb);
   if cr ~=0
      return
   end
   if ~isempty(option)
      param.cons   = zsetfield(param.cons,module_st,option);
   end
elseif ~isempty(strmatch(module,{'coef_all','onestep','first','sources'},'exact'))
   % pas de parametres dans ce cas
elseif isempty(zgetfield(param.fonction,module_st))
   fprintf('undefined function for module %s !\n',module);
   return
elseif isempty(option) & (updateparam ~= 0)
   % mise a jour des parametres
   try
      nb = zgetfield(param.nombre,module);
   catch
      nb = 1;
   end
        [param,nbok] =  locmergeparam(param,module_st,zgetfield(param.fonction,module_st),nb);
elseif ~isempty(option)
   param.cons   = zsetfield(param.cons,module_st,option);
end


% recuperation de 1 temps
if length(data.gene.temps) == 1
   datak = data;
   datakp1 = data;
   param.gene.nbt  = 1;
elseif length(data.gene.temps) == 2
   datak = zget1t(data,1);
   datakp1 = zget1t(data,2);
   param.gene.nbt  = 2;
elseif isempty(temps)
   datak = zget1t(data,param.gene.k);
   indnext = min(find(data.gene.temps > datak.gene.temps));
   if ~isempty(indnext)
      datakp1 = zget1t(data,indnext);
   else
      datakp1 = datak;
      datakp1.gene.temps = datak.gene.temps + max(eps,param.gene.dt);
   end
   param.gene.nbt  = 2;
else
   ind = min(find(data.gene.temps >= temps));
   datak = zget1t(data,ind);
   if ind < length(data.gene.temps)
      datakp1 = zget1t(data,ind+1);
   else
      datakp1 = datak;
      datakp1.gene.temps = datak.gene.temps + max(eps,param.gene.dt);
   end
   param.gene.nbt  = 2;
end

% mise a jour de param
if isempty(strmatch(module,{'asser','onestep','first'},'exact'))
   param.gene.k    = 1;
   param.gene.kmin = 1;
   param.gene.kmax = param.gene.nbt;
end

% passage des parametre (scalaire) machines dependant a tous les modules externes
nomm = fieldnames(param.fonction);
indm = strmatch('machine',nomm,'exact');
% compatibilite ascendante
if ~isempty(indm)
   nomm(indm) = [];
   machine_param  = param.cons.machine;
   if ~isempty(machine_param)
      for km = 1:length(nomm)
         consm = zgetfield(param.cons,nomm{km});
         consm = zsetfield(consm,'machine_param',machine_param);
         param.cons  = zsetfield(param.cons,nomm{km},consm);
      end
   end
end
% ajout de la variable standart machine_nom
nomm = fieldnames(param.fonction);
for km = 1:length(nomm)
    consm = zgetfield(param.cons,nomm{km});
    consm = zsetfield(consm,'machine_nom',upper(param.from.machine));
    param.cons  = zsetfield(param.cons,nomm{km},consm);
end

% passage des consigne d'asservissement au modules externes
% (permet d'asservir de grandeurs comme la position des antennes)
nomm = fieldnames(param.fonction);
for km = 1:length(nomm)
    consm = zgetfield(param.cons,nomm{km});
    consm = zsetfield(consm,'machine_cons',datak.cons.asser);
    param.cons  = zsetfield(param.cons,nomm{km},consm);
end

% passage des consigne d'asservissement au modules externes
% (permet d'asservir de grandeurs comme la position des antennes)
stemps.tempsini = datak.gene.temps;     % temps initial de l'intervalle de temps interne en cours
if param.gene.nbt > 1
   stemps.tempsend = datakp1.gene.temps;   % temps final de l'intervalle de temps interne en cours
else
   stemps.tempsend = datak.gene.temps;   % temps final de l'intervalle de temps interne en cours
end
stemps.dtdiff   = param.gene.dt;           % pas de temps du solveur de diffusion
stemps.dtequi   = param.gene.dt;                   % duree ecoule depuyis le calcul du dernier equilibre
stemps.dtasser  = param.gene.dt;                    % pas de temps ineterne pour le calcul des asservissement (peu etre plus grand que dtdiff, selon decoupe)
nomm = fieldnames(param.fonction);
for km = 1:length(nomm)
      consm = zgetfield(param.cons,nomm{km});
      consm = zsetfield(consm,'temps_cons',stemps);
      param.cons  = zsetfield(param.cons,nomm{km},consm);
end

% par default ce n'est pas le premier temps
datak.mode.premiertemps   = premiertemps;
datakp1.mode.premiertemps = premiertemps;


% les modules doivent ecrire les infos
param.gene.verbose  =1;
setappdata(0,'zverbose',1);

% pour l'appel des sources
proto          = zsourceproto(param.gene.nbrho);

% changement des chemins de sortie
memgene                 = param.gene;
home                    = getenv('HOME');
if ~isempty(nomf)
   rcfile  = strcat('test_',module,'_',nomf,'@',datestr(now,30));
else
   rcfile  = strcat('test_',module,'@',datestr(now,30));
end
rcfile =fullfile('.',rcfile);
param.gene.rcfile = rcfile;

if isdir(fullfile(home,'zineb'))
   if isdir(fullfile(home,'zineb','data'))
      param.gene.origine = fullfile(home,'zineb','data',rcfile);
      param.gene.file         = strcat(param.gene.origine,'_resultat');
      if isdir(fullfile(home,'zineb','data','rapsauve'))
         param.gene.rapsauve    = fullfile(home,'zineb','data','rapsauve',rcfile);
      else
         param.gene.rapsauve    = strcat(param.gene.origine,'_resultat_rapsauve');
      end
   else
      param.gene.origine   = rcfile;
      param.gene.file         = strcat(param.gene.origine,'_resultat');
      param.gene.rapaauve     = strcat(param.gene.origine,'_resultat_rapsauve');
   end
else
   param.gene.origine      = rcfile;
   param.gene.file         = strcat(param.gene.origine,'_resultat');
   param.gene.rpasauve     = strcat(param.gene.origine,'_resultat_rapsauve');
end

% definition du repertoire de travail
setappdata(0,'CRONOS_WORK_DIR',fileparts(param.gene.file));

% appel du module
try
   if debug_onoff
      dbstop if all error
   elseif ~isempty(dbstatus)
      disp('debug mode is :')
      dbstatus
   end

   switch module
   case 'impur'
      fprintf('Call of %s (impurity) with parameters :\n',param.fonction.impur);
      parameters = param.cons.impur
      cons = param.cons.impur;
      cons.premiertemps = premiertemps;
      % passage du flag d'initilisation
      [impur,memoire]  = ...
      feval(param.fonction.impur,cons, datak.geo,datak.equi,datak.bord, ...
            datak.cons.zeffm,datak.cons.nhnd,datak.prof,datak.neo, ...
            param.phys,param.compo,param.gene, ...
            param.memoire.impur,datak.source.n,datak.impur, ....
            param.fonction.neo,param.cons.neo,datak.coef,param.memoire.neo);
      dbclear all
      memoire.t = param.gene.t;
      fprintf('difference in impur structure :\n');
      zcompstruct(impur,datakp1.impur,tolerance);
      fprintf('difference in memoire.impur structure :\n');
      zcompstruct(memoire,memin.impur,tolerance);
      if plotonoff
         zplotstruct(impur,datakp1.impur,'data_impur',param.fonction.impur);
         zplotstruct(memoire,memin.impur,'memoire_impur',param.fonction.impur);
      end
      % recopie pour les sorties
      datakp1.impur = impur;
      param.memoire.impur = memoire;

   case 'equi'
      fprintf('Call of %s (equilibrium) with parameters :\n',param.fonction.equi);
      parameters = param.cons.equi
            if datak.mode.premiertemps == 1
               ip = datak.cons.ip;
            elseif datak.mode.cons.psi  == 0
                  ip = datak.cons.ip;
            else
                  ip=[];
            end
            [equi,mhd_cd,memoire]  =  ...
             feval(param.fonction.equi,param.memoire.equi,param.cons.equi,datak.cons.asser.pfcur, ...
                              datak.geo,datak.equi,datak.prof, ...
               param.phys,ip,param.gene.x,1);
      dbclear all
      param.memoire.equi = memoire;      
      if isfield(memin.equi,'data')
      		memin.equi.data.equifree_last = [];
       		memin.equi.data.equifree.param.machine_param = [];
      end
      if isfield(memoire,'data')
      		memoire.data.equifree_last = [];
      		memoire.data.equifree.param.machine_param = [];
      end
     
      fprintf('difference in equi structure :\n');
      equi.errit = datakp1.equi.errit;
      zcompstruct(equi,datakp1.equi,tolerance);
      fprintf('difference in prof.mhd_cd structure :\n');
      zcompstruct(mhd_cd,datakp1.prof.mhd_cd,tolerance);
       fprintf('difference in param.memoire.equi structure :\n');
      zcompstruct(memoire,memin.equi,tolerance);
      if plotonoff
         zplotstruct(equi,datakp1.equi,'data_equi',param.fonction.equi);
         zplotstruct(mhd_cd,datakp1.prof.mhd_cd,'data_prof_mhd_cd',param.fonction.equi);
         zplotstruct(memoire,memin.equi,'param_memoire_equi',param.fonction.equi);
      end
      % recopie pour les sorties
      datakp1.equi = equi;
      datakp1.prof.mhd_cd  = mhd_cd;
 
   case 'neo'
      fprintf('Call of %s (neoclassical) with parameters :\n',param.fonction.neo);
      parameters = param.cons.neo
      datak.mode.neo = 2;
      if isfield(param.memoire.impur,'etatcharge')
            [cr,neo,memoire] = zneoclass(param.fonction.neo,param.cons.neo, ...
                              param.gene,param.phys,param.compo,datak, ...
                     datak.neo,0,param.memoire.impur.etatcharge,param.memoire.neo);
      else
            [cr,neo,memoire] = zneoclass(param.fonction.neo,param.cons.neo, ...
                              param.gene,param.phys,param.compo,datak, ...
                     datak.neo,0,[],param.memoire.neo);
      end
      dbclear all
      memoire.t = param.gene.t;
      memin.neo.t = param.gene.t;
      fprintf('difference in neo structure :\n');
      zcompstruct(neo,datakp1.neo,tolerance);
      fprintf('difference in memoire.neo structure :\n');
      memoire.data.neodir = '';
      memin.neo.data.neodir ='';
      zcompstruct(memoire,memin.neo,tolerance);
      if plotonoff
         zplotstruct(neo,datakp1.neo,'data_neo',param.fonction.neo);
         zplotstruct(memoire,memin.neo,'memoire_neo',param.fonction.neo);
      end
      % recopie pour les sorties
      datakp1.neo = neo;
      param.memoire.neo = memoire;

   case 'rip'
      fprintf('Call of %s (ripple) with parameters :\n',param.fonction.rip);
      parameters = param.cons.rip
      [rip,nrip] = feval(param.fonction.rip,param.cons.rip,proto, ...
                   datak.geo,datak.equi,datak.prof,datak.neo, ...
                             datak.impur,datak.bord,param.phys,param.compo,param.gene, ...
              datak.mode.premiertemps);
      dbclear all
      fprintf('difference in rip structure :\n');
      zcompstruct(rip,datakp1.source.rip,tolerance);
      fprintf('difference in source.n.rip structure :\n');
      zcompstruct(nrip,datakp1.source.n.rip,tolerance);
      if plotonoff
         zplotstruct(rip,datakp1.source.rip,'data_source_rip',param.fonction.rip);
         zplotstruct(nrip,datakp1.source.n.rip,'data_source_n_rip',param.fonction.rip);
      end

      % recopie pour les sorties
      datakp1.source.rip = rip;
      datakp1.source.n.rip = nrip;

   case 'mhd_dds'

      fprintf('Call of %s (MHD sawtooth) with parameters :\n',param.fonction.mhd.dds);
      parameters = param.cons.mhd.dds
      [cr,dataknew,memoire] = feval(param.fonction.mhd.dds,param.cons.mhd.dds,datak,param.gene,param.compo, ...
                      param.phys,datak.gene,param.fonction.equi,param.cons.equi,param.memoire.equi,datak.cons.asser.pfcur);

      dbclear all
      fprintf('difference in datak structure after mhd_dds :\n');
      dataknew.gene.temps = datakp1.gene.temps;
      zcompstruct(dataknew,datakp1,tolerance);
      fprintf('difference in memoire structure after mhd_dds :\n');
      zcompstruct(memoire,memin.equi,tolerance);
      datakp1 = dataknew;

   case 'mhd_limite'
      fprintf('Call of %s (MHD trigger) with parameters :\n',param.fonction.mhd.limite);
      parameters = param.cons.mhd.limite
      [cr,evx] = feval(param.fonction.mhd.limite,param.cons.mhd.limite, ...
                  datak,datak,param.gene,param.compo,param.phys);
      dbclear all
      fprintf('difference in evx structure after mhd_limite :\n');
      zcompstruct(evx,datakp1.evx,tolerance);
      datakp1.evx = evx;

   case 'mhd_elm'
      fprintf('Call of %s (MHD Elms) with parameters :\n',param.fonction.mhd.elm);
      parameters = param.cons.mhd.elm
      [cr,dataknew,memoire] = feval(param.fonction.mhd.elm,param.cons.mhd.elm,datak,param.gene,param.compo, ...
                      param.phys,datak.gene,param.fonction.equi,param.cons.equi,param.memoire.equi);
      dbclear all
      fprintf('difference in datak structure after mhd_elm :\n');
      zcompstruct(dataknew,datakp1,tolerance);
      fprintf('difference in memoire structure after mhd_elm :\n');
      zcompstruct(memoire,memin.equi,tolerance);
      datakp1 = dataknew;

   case 'mhd_stab'
      datak.mode.mhd.stab =  2
      fprintf('Call of %s (MHD stability) with parameters :\n',param.fonction.mhd.stab);
      parameters = param.cons.mhd.stab
      [cr,mhd,mhd_cd] = zstabmhd(param.fonction.mhd.stab,param.cons.mhd.stab,param.gene, ...
                                 param.phys,param.compo,datak,datak.equi.mhd);
      dbclear all
      mhd_aux=datakp1.equi.mhd;
      mhd_aux.gamma1=mhd.gamma1;
      mhd_aux.gamma2=mhd.gamma2;
      mhd_aux.gamma3=mhd.gamma3;
      fprintf('difference in equi.mhd structure after mhd_stab (only gamma is significant):\n');
      zcompstruct(mhd_aux,datakp1.equi.mhd,tolerance);
      %fprintf('difference in datak.prof.mhd_cd structure after mhd_stab :\n');
      %zcompstruct(mhd_cd,datakp1.prof.mhd_cd,tolerance);
      if plotonoff
         zplotstruct(mhd,datakp1.equi.mhd,'data_equi_mhd',param.fonction.mhd.stab);
         %zplotstruct(mhd_cd,datakp1.prof.mhd_cd,'data_prof_mhd_cd',param.fonction.mhd.stab);
      end
      % recopie pour les sorties
      datakp1.equi.mhd     = mhd;
      datakp1.prof.mhd_cd  = mhd_cd;

   case 'sources'
      fprintf('Call of all sources :\n');
      data.mode.fci = 2;
      data.mode.fce = 2;
      data.mode.hyb = 2;
      data.mode.idn = 2;
      data.mode.fus = 2;
      data.mode.n0  = 2;
      data.mode.cyclo  = 2;
      [cr,dataknew,memoire1]   = zsources(datak,param.gene,param.compo,param.fonction,param.memoire,param.cons, ...
                                      param.phys,datak.source,param.batch,0);
      dbclear all
      dataknew.gene.temps = datakp1.gene.temps;
      %
      % gestion du temps des memoires
      %
      noms = fieldnames(param.memoire);
      for k=1:length(noms)
         if isfield(param.memoire.(noms{k}),'t')
            param.memoire.(noms{k}).t = param.gene.t;
         end
      end
      memoirebis=memoire1;
      % ajout d'une exception pour pion concernant fr
      reldiff=0.0;
      if isfield(memoire1.fci.data,'output')
        if isfield(memoire1.fci.data.output,'fr')
                valmax=max(max(memin.fci.data.output.fr));
                val=(abs(memoire1.fci.data.output.fr - memin.fci.data.output.fr) )./valmax;
                relaux=max(max(val));
                fprintf(' real value differ by %g of norm: (memoire.data.output.fr)\n',relaux);
                reldiff=max(reldiff,relaux);
                memoirebis.fci.data.output.fr=memin.fci.data.output.fr;
        end
        if isfield(memoire1.fci.data.output,'fl')
                valmax=max(max(memin.fci.data.output.fl));
                val=(abs(memoire1.fci.data.output.fl - memin.fci.data.output.fl) )./valmax;
                relaux=max(max(val));
                fprintf(' real value differ by %g of norm: (memoire.data.output.fl)\n',relaux);
                reldiff=max(reldiff,relaux);
                memoirebis.fci.data.output.fl=memin.fci.data.output.fl;
        end
        if isfield(memoire1.fci.data.output,'fc')
                valmax=max(max(memin.fci.data.output.fc));
                val=(abs(memoire1.fci.data.output.fc - memin.fci.data.output.fc) )./valmax;
                relaux=max(max(val));
                fprintf(' real value differ by %g of norm: (memoire.data.output.fc)\n',relaux);
                reldiff=max(reldiff,relaux);
                memoirebis.fci.data.output.fc=memin.fci.data.output.fc;
        end
        if reldiff < 1.0e-1
                        if reldiff == 0.0
                            fprintf('CRONOSTEST: OK \n');
                        else
                            fprintf('CRONOSTEST: WARNING = %g of norm \n',reldiff);
                        end
        else
                        fprintf('CRONOSTEST: ERROR = %g of norm \n',reldiff);
        end
      end
     
      fprintf('difference in datak structure after sources :\n');
      zcompstruct(dataknew,datakp1,tolerance);
      fprintf('difference in memoire structure after sources :\n');
      zcompstruct(memoirebis,memin,tolerance);
      datakp1 = dataknew;
      param.memoire = memoire1;

   case 'fci'
      fprintf('Call of %s (fci) with parameters :\n',param.fonction.fci);
      parameters = param.cons.fci
            [fci,memoire]  =  feval(param.fonction.fci,param.cons.fci,proto,datak.cons.fci, ...
                                        datak.geo,datak.equi,datak.cons.c,datak.prof,datak.neo, ...
                                        datak.impur,param.phys,param.compo,param.gene,param.memoire.fci);
      dbclear all
      memoire.t = param.gene.t;
      fprintf('difference in fci structure :\n');
      zcompstruct(fci,datakp1.source.fci,tolerance);
      fprintf('difference in memoire.fci structure :\n');
      memoirebis=memoire;
      % ajout d'une exception pour pion concernant fr
      reldiff=0.0;
      if isfield(memoire.data,'output')
        if isfield(memoire.data.output,'fr')
                valmax=max(max(memin.fci.data.output.fr));
                val=(abs(memoire.data.output.fr - memin.fci.data.output.fr) )./valmax;
                relaux=max(max(val));
                fprintf(' real value differ by %g of norm: (memoire.data.output.fr)\n',relaux);
                reldiff=max(reldiff,relaux);
                memoirebis.data.output.fr=memin.fci.data.output.fr;
        end
        if isfield(memoire.data.output,'fl')
                valmax=max(max(memin.fci.data.output.fl));
                val=(abs(memoire.data.output.fl - memin.fci.data.output.fl) )./valmax;
                relaux=max(max(val));
                fprintf(' real value differ by %g of norm: (memoire.data.output.fl)\n',relaux);
                reldiff=max(reldiff,relaux);
                memoirebis.data.output.fl=memin.fci.data.output.fl;
        end
        if isfield(memoire.data.output,'fc')
                valmax=max(max(memin.fci.data.output.fc));
                val=(abs(memoire.data.output.fc - memin.fci.data.output.fc) )./valmax;
                relaux=max(max(val));
                fprintf(' real value differ by %g of norm: (memoire.data.output.fc)\n',relaux);
                reldiff=max(reldiff,relaux);
                memoirebis.data.output.fc=memin.fci.data.output.fc;
        end
        if reldiff < 1.0e-1
                        if reldiff == 0.0
                            fprintf('CRONOSTEST: OK \n');
                        else
                            fprintf('CRONOSTEST: WARNING = %g of norm \n',reldiff);
                        end
        else
                        fprintf('CRONOSTEST: ERROR = %g of norm \n',reldiff);
        end
      end
      zcompstruct(memoirebis,memin.fci,tolerance);
      if plotonoff
         zplotstruct(fci,datakp1.source.fci,'data_source_fci',param.fonction.fci);
         zplotstruct(memoire,memin.fci,'memoire_fci',param.fonction.fci);
      end
      % recopie pour les sorties
      datakp1.source.fci = fci;
      param.memoire.fci = memoire;

   case 'fce'
      if ~isempty(strmatch(lower(param.fonction.fce),{'zecrhinlh'},'exact'))
         fprintf('Call of %s with %s (ECRH in HYB) : no computation with fce key\n',param.fonction.fce,param.fonction.hyb);
         datak1.source.fce   = proto;

      else
         fprintf('Call of %s (fce) with parameters :\n',param.fonction.fce);
         parameters = param.cons.fce
         parameters.luke_mode = 'ecrh';
         [fce,memoire]  =  feval(param.fonction.fce,parameters,proto,datak.cons.fce, ...
                  datak.geo,datak.equi,datak.cons.c,datak.prof,datak.neo, ...
                  datak.impur,param.phys,param.compo,param.gene,param.memoire.fce);
         dbclear all
         memoire.t = param.gene.t;
         memin.fce.t = param.gene.t;
         fprintf('difference in fce structure :\n');
         zcompstruct(fce,datakp1.source.fce,tolerance);
         fprintf('difference in memoire.fce structure :\n');
         zcompstruct(memoire,memin.fce,tolerance);
         if plotonoff
            zplotstruct(fce,datakp1.source.fce,'data_source_fce',param.fonction.fce);
            zplotstruct(memoire,memin.fce,'memoire_fce',param.fonction.fce);
         end
         % recopie pour les sorties
         datakp1.source.fce = fce;
         param.memoire.fce = memoire;
      end

   case 'hyb'
      if ~isempty(strmatch(lower(param.fonction.fce),{'zecrhinlh'},'exact'))
         fprintf('Call of %s (fce + hyb) with parameters :\n',parammixte);
         parameters = param.cons.hyb;
         parameters.lh_param   = param.cons.hyb;
         parameters.ecrh_param = param.cons.fce;
         parameters.luke_mode = 'mixte';
         parameters.lh_param.luke_mode = 'mixte';
         parameters.ecrh_param.luke_mode = 'mixte';
         [hyb,memoire]  =  feval(param.fonction.hyb,parameters,proto, ...
                  datak.cons.hyb,datak.exp.hybspec, ...
                  datak.geo,datak.equi,datak.cons.c,datak.prof,datak.neo, ...
                  datak.impur,param.phys,param.compo,param.gene,param.memoire.hyb);
         dbclear all
         memoire.t = param.gene.t;
         fprintf('difference in hyb structure :\n');
         zcompstruct(hyb,datakp1.source.hyb,tolerance);
         fprintf('difference in memoire.hyb structure :\n');
         zcompstruct(memoire,memin.hyb,tolerance);
         if plotonoff
            zplotstruct(hyb,datakp1.source.hyb,'data_source_hyb',param.fonction.hyb);
            zplotstruct(memoire,memin.hyb,'memoire_hyb',param.fonction.hyb);
         end
         % recopie pour les sorties
         datakp1.source.hyb = hyb;
         param.memoire.hyb = memoire;
      else
         fprintf('Call of %s (hyb) with parameters :\n',param.fonction.hyb);
         parameters = param.cons.hyb
         parameters.luke_mode = 'lh';
         [hyb,memoire]  =  feval(param.fonction.hyb,parameters,proto, ...
                  datak.cons.hyb,datak.exp.hybspec, ...
                  datak.geo,datak.equi,datak.cons.c,datak.prof,datak.neo, ...
                  datak.impur,param.phys,param.compo,param.gene,param.memoire.hyb);
         dbclear all
         memoire.t = param.gene.t;
         fprintf('difference in hyb structure :\n');
         zcompstruct(hyb,datakp1.source.hyb,tolerance);
         fprintf('difference in memoire.hyb structure :\n');
         zcompstruct(memoire,memin.hyb,tolerance);
         if plotonoff
            zplotstruct(hyb,datakp1.source.hyb,'data_source_hyb',param.fonction.hyb);
            zplotstruct(memoire,memin.hyb,'memoire_hyb',param.fonction.hyb);
         end
         % recopie pour les sorties
         datakp1.source.hyb = hyb;
         param.memoire.hyb = memoire;
      end
   case 'idn'
      fprintf('Call of %s (idn) with parameters :\n',param.fonction.idn);
      parameters = param.cons.idn
            [idn,nidn,memoire]  =  feval(param.fonction.idn,param.cons.idn,proto,datak.cons.idn, ...
                                        datak.geo,datak.equi,datak.cons.c,datak.prof,datak.neo, ...
                                        datak.impur,param.phys,param.compo,param.gene,param.memoire.idn);
      dbclear all
      memoire.t = param.gene.t;
      fprintf('difference in idn structure :\n');
      zcompstruct(idn,datakp1.source.idn,tolerance);
      fprintf('difference in source.n.idn structure :\n');
      zcompstruct(nidn,datakp1.source.n.idn,tolerance);
      fprintf('difference in memoire.idn structure :\n');
      zcompstruct(memoire,memin.idn,tolerance);
      if plotonoff
         zplotstruct(idn,datakp1.source.idn,'data_source_idn',param.fonction.idn);
         zplotstruct(memoire,memin.idn,'memoire_idn',param.fonction.idn);
         zplotstruct(nidn,datakp1.source.n.idn,'data_source_n_idn',param.fonction.idn);
      end

      % recopie pour les sorties
      datakp1.source.idn = idn;
      param.memoire.idn = memoire;
      datakp1.source.n.idn = nidn;

   case 'n0'
      fprintf('Call of %s (cold edge neutals) with parameters :\n',param.fonction.n0);
      parameters = param.cons.n0
            [n0,nn0,nn0th,memoire]  = feval(param.fonction.n0,param.cons.n0,proto, ...
                       datak.geo,datak.equi,datak.cons.c, ...
                       datak.prof,datak.neo,datak.impur,datak.bord, ...
                        param.phys,param.compo,param.gene,param.memoire.n0);
      dbclear all
      memoire.t = param.gene.t;
      fprintf('difference in n0 structure :\n');
      zcompstruct(n0,datakp1.source.n0,tolerance);
      fprintf('difference in source.n.n0 structure :\n');
      zcompstruct(nn0,datakp1.source.n.n0,tolerance);
      fprintf('difference in source.n.n0th structure :\n');
      zcompstruct(nn0th,datakp1.source.n.n0th,tolerance);
      fprintf('difference in memoire.n0 structure :\n');
      zcompstruct(memoire,memin.n0,tolerance);
      if plotonoff
         zplotstruct(n0,datakp1.source.n0,'data_source_n0',param.fonction.n0);
         zplotstruct(memoire,memin.n0,'memoire_n0',param.fonction.n0);
         zplotstruct(nn0,datakp1.source.n.n0,'data_source_n_n0',param.fonction.n0);
         zplotstruct(nn0th,datakp1.source.n.n0th,'data_source_n_n0th',param.fonction.n0);
      end
      % recopie pour les sorties
      datakp1.source.n0 = n0;
      param.memoire.n0 = memoire;
      datakp1.source.n.n0 = nn0;
      datakp1.source.n.n0th = nn0th;

   case 'bord'
      fprintf('Call of %s (SOL module) with parameters :\n',param.fonction.bord);
      parameters = param.cons.bord

            pomp = datak.cons.pomp;
            pomp(~isfinite(pomp)) = 0;
            gaz = datak.cons.c;
            gaz(~isfinite(gaz)) =0;

            [bord,memoire]  =  feval(param.fonction.bord,param.cons.bord,datak.bord,pomp,gaz, ...
                                         datak.prof,datak.impur,datak.source, ...
                                         param.compo,param.memoire.bord,datak.equi.vpr,datak.equi.spr, ...
                datak.equi.rhomax,datak.equi.grho2, ...
                                         param.gene.x,param.phys,param.gene.nbg,param.gene.dt);
      dbclear all
      memoire.t = param.gene.t;
      fprintf('difference in bord structure :\n');
      zcompstruct(bord,datakp1.bord,tolerance);
      fprintf('difference in memoire.bord structure :\n');
      zcompstruct(memoire,memin.bord,tolerance);
      if plotonoff
         zplotstruct(bord,datakp1.bord,'data_bord',param.fonction.bord);
         zplotstruct(memoire,memin.bord,'memoire_bord',param.fonction.bord);
      end
      % recopie pour les sorties
      datakp1.bord = bord;
      param.memoire.bord = memoire;

   case 'glacon'
          fprintf('Call of %s (pellets) with parameters :\n',param.fonction.glacon);
          parameters = param.cons.glacon
          datak.mode.glacon = 1;
          datak.cons.glacon = 1;
          [cr,dataknew] = zglacon(datak,param.gene,param.compo,param.fonction,param.cons,param.phys,datak.gene);
          dbclear all
	  dataknew.gene.temps = datakp1.gene.temps;
	  dataknew.gene.neutralite     = datakp1.gene.neutralite;
          fprintf('difference in datak structure after glacon :\n');
          zcompstruct(dataknew,datakp1,tolerance);
          datakp1 = dataknew;

   case 'fus'
      fprintf('Call of %s (fus) with parameters :\n',param.fonction.fus);
      parameters = param.cons.fus
            [fus,nfus,memoire] = feval(param.fonction.fus,param.cons.fus,proto, ...
                        datak.source.n.fus,datak.geo,datak.equi,datak.cons.c, ...
              datak.prof,datak.neo,datak.impur, ...
              param.phys,param.compo,param.gene,param.memoire.fus);
      dbclear all
      memoire.t = param.gene.t;
      fprintf('difference in fus structure :\n');
      zcompstruct(fus,datakp1.source.fus,tolerance);
      fprintf('difference in source.n.fus structure :\n');
      zcompstruct(nfus,datakp1.source.n.fus,tolerance);
      fprintf('difference in memoire.fus structure :\n');
      memoire.data.jobfile='';
      memin.fus.data.jobfile='';
      zcompstruct(memoire,memin.fus,tolerance);
      if plotonoff
         zplotstruct(fus,datakp1.source.fus,'data_source_fus',param.fonction.fus);
         zplotstruct(memoire,memin.fus,'fus',param.fonction.fus);
         zplotstruct(nfus,datakp1.source.n.fus,'data_source_n_fus',param.fonction.fus);
      end
      % recopie pour les sorties
      datakp1.source.fus = fus;
      param.memoire.fus = memoire;
      datakp1.source.n.fus = nfus;

   case 'cyclo'
      fprintf('Call of %s (cyclo) with parameters :\n',param.fonction.cyclo);
      parameters = param.cons.cyclo
      [cyclo,memoire]  = ...
         feval(param.fonction.cyclo,param.cons.cyclo, ...
         datak.geo,datak.equi,datak.prof,datak.neo,datak.impur, ...
         param.phys,param.compo,param.gene,param.memoire.cyclo);
      dbclear all
      memoire.t = param.gene.t;
      fprintf('difference in cyclo structure :\n');
      zcompstruct(cyclo,datakp1.source.cyclo,tolerance);
      fprintf('difference in memoire.cyclo structure :\n');
      %zcompstruct(memoire,memin.cyclo,tolerance);
      if plotonoff
         zplotstruct(cyclo,datakp1.source.cyclo,'data_source_cyclo',param.fonction.cyclo);
         zplotstruct(memoire,memin.cyclo,'memoire_cyclo',param.fonction.cyclo);
      end
      % recopie pour les sorties
      datakp1.source.cyclo = cyclo;
      param.memoire.cyclo = memoire;

   case 'coefa'
      if ~isempty(param.fonction.coefa)
            fprintf('Call of %s (transport coefficient module A) with parameters :\n',param.fonction.coefa);
            parameters = param.cons.coefa
            coef  =  feval(param.fonction.coefa,param.cons.coefa,datak.coef, ...
                           datak.geo,datak.equi,datak.impur,datak.prof,datak.neo, ...
                           param.phys,param.compo,param.gene.nbrho,param.gene.nbg, ...
               param.gene.x,datak.source);
      end
      dbclear all
      fprintf('difference in coef structure after coefa :\n');
      zcompstruct(coef,datakp1.coef,tolerance);
      if plotonoff
         zplotstruct(coef,datakp1.coef,'data_coef_a',param.fonction.coefa);
      end
      % recopie pour les sorties
      datakp1.coef = coef;

   case 'coefb'
      if ~isempty(param.fonction.coefb)
            fprintf('Call of %s (transport coefficient module B) with parameters :\n',param.fonction.coefb);
            parameters = param.cons.coefb
            coef  =  feval(param.fonction.coefb,param.cons.coefb,datak.coef, ...
                           datak.geo,datak.equi,datak.impur,datak.prof,datak.neo, ...
                           param.phys,param.compo,param.gene.nbrho,param.gene.nbg, ...
               param.gene.x,datak.source);
      end
      dbclear all
      fprintf('difference in coef structure after coefb :\n');
      zcompstruct(coef,datakp1.coef,tolerance);
      if plotonoff
         zplotstruct(coef,datakp1.coef,'data_coef_b',param.fonction.coefb);
      end
      % recopie pour les sorties
      datakp1.coef = coef;

   case 'coefc'
      if ~isempty(param.fonction.coefc)
            fprintf('Call of %s (transport coefficient module C) with parameters :\n',param.fonction.coefc);
            parameters = param.cons.coefc
            coef  =  feval(param.fonction.coefc,param.cons.coefc,datak.coef, ...
                           datak.geo,datak.equi,datak.impur,datak.prof,datak.neo, ...
                           param.phys,param.compo,param.gene.nbrho,param.gene.nbg, ...
               param.gene.x,datak.source);
      end
      dbclear all
      fprintf('difference in coef structure after coefc :\n');
      zcompstruct(coef,datakp1.coef,tolerance);
      if plotonoff
         zplotstruct(coef,datakp1.coef,'data_coef_c',param.fonction.coefc);
      end
      % recopie pour les sorties
      datakp1.coef = coef;

   case 'coefd'
      if ~isempty(param.fonction.coefd)
            fprintf('Call of %s (transport coefficient module D) with parameters :\n',param.fonction.coefd);
            parameters = param.cons.coefd
            coef  =  feval(param.fonction.coefd,param.cons.coefd,datak.coef, ...
                           datak.geo,datak.equi,datak.impur,datak.prof,datak.neo, ...
                           param.phys,param.compo,param.gene.nbrho,param.gene.nbg, ...
               param.gene.x,datak.source);
      end
      dbclear all
      fprintf('difference in coef structure after coefd :\n');
      zcompstruct(coef,datakp1.coef,tolerance);
      if plotonoff
         zplotstruct(coef,datakp1.coef,'data_coef_d',param.fonction.coefd);
      end
      % recopie pour les sorties
      datakp1.coef = coef;

   case 'coefe'
      if ~isempty(param.fonction.coefe)
            fprintf('Call of %s (transport coefficient module E) with parameters :\n',param.fonction.coefe);
            parameters = param.cons.coefe
            coef  =  feval(param.fonction.coefe,param.cons.coefe,datak.coef, ...
                           datak.geo,datak.equi,datak.impur,datak.prof,datak.neo, ...
                           param.phys,param.compo,param.gene.nbrho,param.gene.nbg, ...
               param.gene.x,datak.source);
      end
      dbclear all
      fprintf('difference in coef structure after coefe :\n');
      zcompstruct(coef,datakp1.coef,tolerance);
      if plotonoff
         zplotstruct(coef,datakp1.coef,'data_coef_e',param.fonction.coefe);
      end
      % recopie pour les sorties
      datakp1.coef = coef;

   case 'coeff'
      if ~isempty(param.fonction.coeff)
            fprintf('Call of %s (transport coefficient module F) with parameters :\n',param.fonction.coeff);
            parameters = param.cons.coeff
            coef  =  feval(param.fonction.coeff,param.cons.coeff,datak.coef, ...
                           datak.geo,datak.equi,datak.impur,datak.prof,datak.neo, ...
                           param.phys,param.compo,param.gene.nbrho,param.gene.nbg, ...
               param.gene.x,datak.source);
      end
      dbclear all
      fprintf('difference in coef structure after coeff :\n');
      zcompstruct(coef,datakp1.coef,tolerance);
      if plotonoff
         zplotstruct(coef,datakp1.coef,'data_coef_f',param.fonction.coeff);
      end
      % recopie pour les sorties
      datakp1.coef = coef;

   case 'coef_all'
      fprintf('Call of  all transport coefficient modules  :\n');
         if param.gene.modecoef  ==0
            [cr,dataknew]   = zcoefficient(datak,param.gene,param.compo,param.fonction,param.cons,param.phys);
         else
            [cr,dataknew]   = zcoefficient_conv(datak,param.gene,param.compo,param.fonction,param.cons,param.phys);
         end
            dbclear all
            fprintf('difference in datak.coef structure after coef_all :\n');
            zcompstruct(dataknew.coef,datakp1.coef,tolerance);
            datakp1 = dataknew;

   case 'plot'
      disp('not yet implemented')

   case 'asser'
      fprintf('Call of %s (feedback control module) with parameters :\n',param.fonction.asser);
      % mise a zero des derivees temporelles
      datak.equi.drhomaxdt      = zeros(size(datak.equi.drhomaxdt));
      datak.equi.dvprdt         = zeros(size(datak.equi.dvprdt));
      datak.equi.dsprdt         = zeros(size(datak.equi.dsprdt));
      datak.equi.dphidt         = zeros(size(datak.equi.dphidt));
      datak.prof.dpsidt         = zeros(size(datak.prof.dpsidt));
      datak.prof.dnedt          = zeros(size(datak.prof.dpsidt));
      datak.prof.dpedt          = zeros(size(datak.prof.dpsidt));
      datak.prof.dpiondt        = zeros(size(datak.prof.dpsidt));
      datak.prof.daedt          = zeros(size(datak.prof.dpsidt));
      datak.prof.dflucedt       = zeros(size(datak.prof.dpsidt));
      datak.prof.dfluciondt     = zeros(size(datak.prof.dpsidt));
      datak.prof.drotdt         = zeros(size(datak.prof.dpsidt));
      datak.prof.dpsidt3p       = zeros(size(datak.prof.dpsidt));
      datak.prof.dnedt3p        = zeros(size(datak.prof.dpsidt));
      datak.prof.dpedt3p        = zeros(size(datak.prof.dpsidt));
      datak.prof.dpiondt3p      = zeros(size(datak.prof.dpsidt));
      datak.prof.daedt3p        = zeros(size(datak.prof.dpsidt));
      datak.prof.dflucedt3p     = zeros(size(datak.prof.dpsidt));
      datak.prof.dfluciondt3p   = zeros(size(datak.prof.dpsidt));
      datak.prof.drotdt3p       = zeros(size(datak.prof.dpsidt));
      parameters = param.cons.asser
      datak.mode.asser = 2
      dataplus         = datak;
      dataplus.gene.temps = param.gene.dt + dataplus.gene.temps;
            [cr,consasser,asser] = zasser(param.fonction.asser,param.cons.asser,param.gene,param.compo, ...
                                          param.phys,datak,datak.cons,param.asser,param.gene.dt);
      dbclear all
      fprintf('difference in asser structure after feedback call:\n');
      zcompstruct(asser,assermem,tolerance);
      fprintf('difference in data.cons structure after feedback call:\n');
      zcompstruct(consasser,datakp1.cons,tolerance);
      if plotonoff
         zplotstruct(asser,assermem,'param_asser',param.fonction.asser);
         zplotstruct(consasser,datakp1.cons,'data_cons',param.fonction.asser);
      end
           % recopie des sorties
         datakp1.cons    = consasser;
      param.assermem  = asser;

   case 'machine'
      disp('this is not a true module');

   case 'post'
      disp('not yet implemented')

   case 'onestep'
      fprintf('Call of Cronos for one time step :\n');
      if ~isempty(datakp1)
         % mise a zero des derivees temporelles
         datak.equi.drhomaxdt      = zeros(size(datak.equi.drhomaxdt));
         datak.equi.dvprdt         = zeros(size(datak.equi.dvprdt));
         datak.equi.dsprdt         = zeros(size(datak.equi.dsprdt));
         datak.equi.dphidt         = zeros(size(datak.equi.dphidt));
         datak.prof.dpsidt         = zeros(size(datak.prof.dpsidt));
         datak.prof.dnedt          = zeros(size(datak.prof.dpsidt));
         datak.prof.dpedt          = zeros(size(datak.prof.dpsidt));
         datak.prof.dpiondt        = zeros(size(datak.prof.dpsidt));
         datak.prof.daedt          = zeros(size(datak.prof.dpsidt));
         datak.prof.dflucedt       = zeros(size(datak.prof.dpsidt));
         datak.prof.dfluciondt     = zeros(size(datak.prof.dpsidt));
         datak.prof.drotdt         = zeros(size(datak.prof.dpsidt));
         datak.prof.dpsidt3p       = zeros(size(datak.prof.dpsidt));
         datak.prof.dnedt3p        = zeros(size(datak.prof.dpsidt));
         datak.prof.dpedt3p        = zeros(size(datak.prof.dpsidt));
         datak.prof.dpiondt3p      = zeros(size(datak.prof.dpsidt));
         datak.prof.daedt3p        = zeros(size(datak.prof.dpsidt));
         datak.prof.dflucedt3p     = zeros(size(datak.prof.dpsidt));
         datak.prof.dfluciondt3p   = zeros(size(datak.prof.dpsidt));
         datak.prof.drotdt3p       = zeros(size(datak.prof.dpsidt));

         % speudo jeu de donnees
         dte         = sqrt(param.split.dtmin .* param.split.dtmax);
         param.gene.t       = datak.gene.temps;
         datakn             = datak;
         if ~isempty(timestep)
            datakn.gene.temps  = datak.gene.temps + timestep;
         else
            datakn.gene.temps  = datakp1.gene.temps;
         end
         param.gene.dt      = max(eps,datakn.gene.temps - datak.gene.temps);
         % passage des consigne d'asservissement au modules externes
         % (permet d'asservir de grandeurs comme la position des antennes)
         nomm = fieldnames(param.fonction);
         for km = 1:length(nomm)
            consm = zgetfield(param.cons,nomm{km});
            consm = zsetfield(consm,'machine_cons',datakp1.cons.asser);
            param.cons  = zsetfield(param.cons,nomm{km},consm);
         end
	
	 % poure le calcul
	 datakn.mode.premiertemps   = premiertemps;

         % appel du module decoupe dans le temps
         [cr,evenement,dtc,dataknew,asser,dte,memoire] = zcoupetemps(datak,datakn,param,dte);
         dbclear all
         % les donnees suivante doivent etre comparee a leur champs de reference; on les mets a 0
         %dataknew.gene.neutralite = 0;
         dataknew.gene.neutralite     = dataknew.gene.neutralite ./ max(1,dataknew.gene.netot);
         if abs(dataknew.gene.neutralite) < 1.e-9 
             dataknew.gene.neutralite = 0.0
         end
         dataknew.equi.drhomaxdt      = zeros(size(datak.equi.drhomaxdt));
         dataknew.equi.dvprdt         = zeros(size(datak.equi.dvprdt));
         dataknew.equi.dsprdt         = zeros(size(datak.equi.dsprdt));
         dataknew.equi.dphidt         = zeros(size(datak.equi.dphidt));
         dataknew.prof.dpsidt         = zeros(size(datak.prof.dpsidt));
         dataknew.prof.dnedt          = zeros(size(datak.prof.dpsidt));
         dataknew.prof.dpedt          = zeros(size(datak.prof.dpsidt));
         dataknew.prof.dpiondt        = zeros(size(datak.prof.dpsidt));
         dataknew.prof.daedt          = zeros(size(datak.prof.dpsidt));
         dataknew.prof.dflucedt       = zeros(size(datak.prof.dpsidt));
         dataknew.prof.dfluciondt     = zeros(size(datak.prof.dpsidt));
         dataknew.prof.drotdt         = zeros(size(datak.prof.dpsidt));
         dataknew.prof.dpsidt3p       = zeros(size(datak.prof.dpsidt));
         dataknew.prof.dnedt3p        = zeros(size(datak.prof.dpsidt));
         dataknew.prof.dpedt3p        = zeros(size(datak.prof.dpsidt));
         dataknew.prof.dpiondt3p      = zeros(size(datak.prof.dpsidt));
         dataknew.prof.daedt3p        = zeros(size(datak.prof.dpsidt));
         dataknew.prof.dflucedt3p     = zeros(size(datak.prof.dpsidt));
         dataknew.prof.dfluciondt3p   = zeros(size(datak.prof.dpsidt));
         dataknew.prof.drotdt3p       = zeros(size(datak.prof.dpsidt));

         %
         % gestion du temps des memoires
         %
         noms = fieldnames(param.memoire);
         for k=1:length(noms)
            if isfield(param.memoire.(noms{k}),'t')
               param.memoire.(noms{k}).t = param.gene.t;
            end
         end
 
        reldiff=0.0;
        memoirebis=memoire;
        if isfield(memoire.fci.data,'output')
           fprintf('difference in memoire.fci structure :\n');
           if isfield(memoire.fci.data.output,'fr')
                   valmax=max(max(memin.fci.data.output.fr));
                   val=(abs(memoire.fci.data.output.fr - memin.fci.data.output.fr) )./valmax;
                   relaux=max(max(val));
                   fprintf(' real value differ by %g of norm: (memoire.data.output.fr)\n',relaux);
                   reldiff=max(reldiff,relaux);
                   memoirebis.fci.data.output.fr=memin.fci.data.output.fr;
           end
           if isfield(memoire.fci.data.output,'fl')
                   valmax=max(max(memin.fci.data.output.fl));
                   val=(abs(memoire.fci.data.output.fl - memin.fci.data.output.fl) )./valmax;
                   relaux=max(max(val));
                   fprintf(' real value differ by %g of norm: (memoire.data.output.fl)\n',relaux);
                   reldiff=max(reldiff,relaux);
                   memoirebis.fci.data.output.fl=memin.fci.data.output.fl;
           end
           if isfield(memoire.fci.data.output,'fc')
                   valmax=max(max(memin.fci.data.output.fc));
                   val=(abs(memoire.fci.data.output.fc - memin.fci.data.output.fc) )./valmax;
                   relaux=max(max(val));
                   fprintf(' real value differ by %g of norm: (memoire.data.output.fc)\n',relaux);
                   reldiff=max(reldiff,relaux);
                   memoirebis.data.output.fc=memin.fci.data.output.fc;
           end
           if reldiff < 1.0e-1
                  if reldiff == 0.0
                         fprintf('CRONOSTEST: OK \n');
                  else
                         fprintf('CRONOSTEST: WARNING = %g of norm \n',reldiff);
                  end
           else
                  fprintf('CRONOSTEST: ERROR = %g of norm \n',reldiff);
           end
        end
	if isfield(memin.neo.data,'neodir')
		memoirebis.neo.data.neodir = memin.neo.data.neodir;
	end 

 

               fprintf('difference in datak structure after one time step :\n');
               zcompstruct(dataknew,datakp1,tolerance);
               fprintf('difference in asser structure after one time step :\n');
               zcompstruct(asser,assermem,tolerance);
               fprintf('difference in memoire structure after one time step :\n');
               zcompstruct(memoirebis,memin,tolerance);

         % copie
         datakp1 = dataknew;
         % memorisation de la derniere valeur de la structure des assservissemnts
         param.assermem   = asser;
         % memorisation de la derniere valeur calculee des sources
         param.memoire   = memoire;
      else
         disp('undefined second time data')
      end

   case 'first'
      fprintf('Call of CRONOS initialisation phase :\n');
      data_ini  = datak;
      param_ini = param;
      param_ini.gene.nbt = 1;
      param_ini.gene.k = 1;
      param_ini.gene.kmin = 1;
      param_ini.gene.kmax = 1;
      param_ini.gene.t = datak.gene.temps;
      dtmem = max(1e-3,param.gene.dt);
      datakp1.gene.temps = datak.gene.temps + dtmem;
      [param_ini,data_ini]=zadd1t(param_ini,data_ini,datakp1);
      [cr,data_ini,param_out]=zpremiertemps(data_ini,param_ini);
      data_out = zget1t(data_ini,1);
      data_out.gene.temps = datakp1.gene.temps;
      data_out.gene.datation = datakp1.gene.datation;
      data_out.gene.neutralite = data_out.gene.neutralite ./ max(1,data_out.gene.netot);
      if abs(data_out.gene.neutralite) < 1.e-9
           data_out.gene.neutralite = 0.0
      end

      param_ini.memoire = param_out.memoire;
      param_ini.asser   = param_out.asser;
      %
      % gestion du temps des memoires
      %
      noms = fieldnames(param_out.cons);
      for k=1:length(noms)
         if isfield(param_out.cons.(noms{k}),'machine_param') & ...
            isfield(param_ini.cons.(noms{k}),'machine_param')
            param_out.cons.(noms{k}).machine_param = param_ini.cons.(noms{k}).machine_param;
         end
      end
      noms = fieldnames(param_out.cons);
      for k=1:length(noms)
         if isfield(param_out.cons.(noms{k}),'temps_cons') & ...
            isfield(param_ini.cons.(noms{k}),'temps_cons')
            param_out.cons.(noms{k}).temps_cons = param_ini.cons.(noms{k}).temps_cons;
         end
      end
      dbclear all
      reldiff=0.0;
      memoirebis=param_out.memoire.fci;
      memin =param_ini.memoire;
      if isfield(memoirebis.data,'output')
           fprintf('difference in memoire.fci structure :\n');
           if isfield(memoirebis.data.output,'fr')
                   valmax=max(max(memin.fci.data.output.fr));
                   val=(abs(memoirebis.data.output.fr - memin.fci.data.output.fr) )./valmax;
                   relaux=max(max(val));
                   fprintf(' real value differ by %g of norm: (memoire.data.output.fr)\n',relaux);
                   reldiff=max(reldiff,relaux);
                   memoirebis.data.output.fr=memin.fci.data.output.fr;
           end
           if isfield(memoirebis.data.output,'fl')
                   valmax=max(max(memin.fci.data.output.fl));
                   val=(abs(memoirebis.data.output.fl - memin.fci.data.output.fl) )./valmax;
                   relaux=max(max(val));
                   fprintf(' real value differ by %g of norm: (memoire.data.output.fl)\n',relaux);
                   reldiff=max(reldiff,relaux);
                   memoirebis.data.output.fl=memin.fci.data.output.fl;
           end
           if isfield(memoirebis.data.output,'fc')
                   valmax=max(max(memin.fci.data.output.fc));
                   val=(abs(memoirebis.data.output.fc - memin.fci.data.output.fc) )./valmax;
                   relaux=max(max(val));
                   fprintf(' real value differ by %g of norm: (memoire.data.output.fc)\n',relaux);
                   reldiff=max(reldiff,relaux);
                   memoirebis.data.output.fc=memin.fci.data.output.fc;
           end
           if reldiff < 1.0e-1
                  if reldiff == 0.0
                         fprintf('CRONOSTEST: OK \n');
                  else
                         fprintf('CRONOSTEST: WARNING = %g of norm \n',reldiff);
                  end
           else
                  fprintf('CRONOSTEST: ERROR = %g of norm \n',reldiff);
           end
      end
      param_out.memoire.fci=memoirebis; 
      fprintf('difference in datak structure after CRONOS initilisation :\n');
      zcompstruct(data_out,datakp1,tolerance);
      fprintf('difference in param structure after CRONOS initilisation :\n');
      zcompstruct(param_out,param_ini,tolerance);
      datakp1 = data_out;
      param = param_out;

   otherwise
      disp('bad way');
   end

   % memorisation du resultat
   param.gene.nbt = length(datak.gene.temps);
   param.gene.k   = 1;
   param.gene.kmin = 1;
   param.gene.kmax = param.gene.nbt;
   if datakp1.gene.temps(1) == datak.gene.temps(1)
      datakp1.gene.temps(1) = datak.gene.temps(1) + 1;
      param.gene.tmax = datak.gene.temps(end);
   end
   [param,datak]=zadd1t(param,datak,datakp1);
   dbclear all

catch
   if ~isempty(strmatch(module,{'coef_all','onestep','first','sources'},'exact'))
      fprintf('Error during execution of %s (%s) : \n', ...
                module,lasterr);
   else
      if ~isempty(findstr(module,'_'))
         [r1,tt]=strtok(module,'_');
         [r2,tt]=strtok(tt,'_');	 
      	 fprintf('Error during execution of %s (%s) : \n%s\n', ...
               	 param.fonction.(r1).(r2),module,lasterr);
     	     
      else
     	 fprintf('Error during execution of %s (%s) : \n%s\n', ...
               	 getfield(param.fonction,module),module,lasterr);
      end
   end
   dbclear all
   if debug_onoff
      disp('Execution stop in zineb1test =>>');
      keyboard
   end
end

% menage
if clean == 1
   [s,t] = unix(sprintf('rm -rf %s*',param.gene.origine));
   [s,t] = unix(sprintf('rm -rf %s*',param.gene.file));
   [s,t] = unix(sprintf('rm -rf %s*',param.gene.rapsauve));
end

% retour au nom d'origine
param.gene.origine    = memgene.origine;
param.gene.file       = memgene.file;
param.gene.rapsauve   = memgene.rapsauve;

% redefinition du repertoire de travail
setappdata(0,'CRONOS_WORK_DIR',fileparts(param.gene.file));


% recopie des parametres communs
function [param,nbok,cr] =  locmergeparam(param,module,nomf,nb)

% valeur de depart
nbok   = 0;
cr     = 0;

% lecture des donnees
try
   par    = feval(nomf,nb);disp('-----------------')
catch
   fprintf('Module %s is unknown or had made an error :\n%s\n',nomf,lasterr);
end
if isfield(par,'valeur');
   valeur = par.valeur;
   if isempty(valeur)
      return
   end
else
   % fonction sans parametre
   return
end

try
   cons = getfield(param.cons,module);
catch
   cons = [];
end

% si pas de parametres existant
if isempty(cons)
   param.cons = zsetfield(param.cons,module,valeur);
   return
end

% boucle sur les champs
noms = fieldnames(valeur);
for k = 1:length(noms)
   nomc = noms{k};
   % si le champ existe deja
   if isfield(cons,nomc)
      valcons = getfield(cons,nomc);
      valval  = getfield(valeur,nomc);
      % si le type est commun
      if (iscell(valval) & iscell(valcons)) | ...
         (ischar(valval) & ischar(valcons)) | ...
         (isnumeric(valval) & isnumeric(valcons))
         % si memes dimensions
         if all(size(valval) == size(valcons)) & isfield(par.borne,nomc)
            borne  = getfield(par.borne,nomc);
            % si intervalle compatible
            ok = inborne(valcons,borne);
            if any(ok)
               valval(find(ok)) = valcons(find(ok));
               eval(sprintf('valeur.%s = valval;',nomc));
               nbok = nbok + 1;
            end
         end
      end
   end
end

% ecriture dans le workspace
param.cons = zsetfield(param.cons,module,valeur);

function ok = inborne(val,borne)

% code retour
ok = ones(1,length(val));

for k = 1:length(val)
   % selon le type
   if iscell(val)
      vc  = val{k};
      if iscell(borne)
         if ischar(vc)
            if ~any(strmatch(vc,borne,'exact'))
               ok(k) = 0;
            end
         else
            borne = zcell2mat(borne);
            if ~any(vc == borne)
               ok(k) = 0;
            end
         end
      else
         if vc < min(borne)
            ok(k) = 0;
         end
         if vc > max(borne)
            ok(k) = 0;
         end
      end
   elseif isnumeric(val)
      vc  = val(k);
      if iscell(borne)
         borne = zcell2mat(borne);
         if ~any(vc == borne)
            ok(k) = 0;
         end
      else
         if vc < min(borne)
            ok(k) = 0;
         end
         if vc > max(borne)
            ok(k) = 0;
         end
      end
   else
      ok(k) = 0;
   end
end
