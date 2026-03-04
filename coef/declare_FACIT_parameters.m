function param = declare_FACIT_parameters

valeur.ION_REGIME	= 1;
type.ION_REGIME	    	= 'real';
borne.ION_REGIME	= {0,1};
defaut.ION_REGIME	= 1;
info.ION_REGIME		= 'Collisionality model\n 0: PS only \n 1: all collisionality regimes';

valeur.OPT_ASYM	= 0;
type.OPT_ASYM	    	= 'real';
borne.OPT_ASYM	= {0,1};
defaut.OPT_ASYM	= 0;
info.OPT_ASYM		= 'Poloidal asymmetries on/off';

valeur.LAB_GEO	= 0;
type.LAB_GEO   	= 'real';
borne.LAB_GEO	= {0,1};
defaut.LAB_GEO	= 0;
info.LAB_GEO	= 'Geometry model. \n0 : circular model \n 1 : arbitrary geometry (requires filling dpsidx, FV, BV, RV from Grad-Shafranov equilibrium)';

valeur.MULTIS	= 0;
type.MULTIS   	= 'real';
borne.MULTIS	= {0,1,2};
defaut.MULTIS	= 0;
info.MULTIS	= 'Accounts for multi species collisions in FACIT with additional simulations. This is also automatically triggered the zeff contribution of the light impurity is above 0.1';

valeur.PLOT	= 0;
type.PLOT   	= 'real';
borne.PLOT	= {0,1};
defaut.PLOT	= 0;
info.PLOT	= 'Plot FACIT inputs/outputs for debugging purposes';


interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet

param.valeur     = valeur;
param.type       = type;
param.borne      = borne;
param.defaut     = defaut;
param.info       = info;
param.interface  = interface;
param.description = 'GUI for FACIT';   % description (une ligne) de la fonction
param.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
param.gui      ='';                             % nom de l'interface graphique specifique si elle existe
param.controle = '';                        % nom de la fonction de controle des valeurs si elle existe

