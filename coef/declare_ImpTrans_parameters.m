%Function used to define Impurity Transport parameters in METIS
function param = declare_ImpTrans_parameters

% METIS implementation parameters
valeur.auto_cons	= 0;
type.auto_cons		= 'real';
borne.auto_cons		= {0,1};
defaut.auto_cons	= 0;
info.auto_cons		= '0 returns particle fluxes and transport coefficients \n1, returns steady state of transport equation';

valeur.neoclass_transport	= '0';
type.neoclass_transport		= 'string';
borne.neoclass_transport		= {'0','FACIT','NEO'};
defaut.neoclass_transport	= '0';
info.neoclass_transport		= 'if = 0, no neoclassical contribution';

valeur.turb_transport	= '0';
type.turb_transport		= 'string';
borne.turb_transport		= {'0','QLKZ_e','QLKZ_imp'};
defaut.turb_transport	= '0';
info.turb_transport		= 'if = 0, no turbulent contribution \n if =QLKZ_e, only works without neoclassical transport (computes steady state for main electrons and main ions. \n if=QLKZ_imp, computes impurity transport coefficients from QLKZ.';

valeur.inter_nbloop_max	= 200;
type.inter_nbloop_max		= 'real';
borne.inter_nbloop_max	= [0,500];
defaut.inter_nbloop_max	= 200;
info.inter_nbloop_max		= 'Maximum number of interation in the internal loop';

valeur.ext_nbloop_max	= 1;
type.ext_nbloop_max		= 'real';
borne.ext_nbloop_max	= [0,50];
defaut.ext_nbloop_max	= 1;
info.ext_nbloop_max		= 'Maximum number of interation in the external loop';


interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet

param.valeur     = valeur;
param.type       = type;
param.borne      = borne;
param.defaut     = defaut;
param.info       = info;
param.interface  = interface;
param.description = 'GUI for ImpTrans';   % description (une ligne) de la fonction
param.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
param.gui      ='';                             % nom de l'interface graphique specifique si elle existe
param.controle = '';                        % nom de la fonction de controle des valeurs si elle existe

