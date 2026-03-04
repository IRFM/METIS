function param = declare_NEO_parameters

valeur.SIM_MODEL	= 1;
type.SIM_MODEL	    	= 'real';
borne.SIM_MODEL		= {0,1,2,3,4};
defaut.SIM_MODEL	= 1;
info.SIM_MODEL		= 'Parameter which selects whether to determine the neoclassical transport from analytic theory or from numerical solution of the DKE.\nSIM_MODEL = 0: analytic theory only \nSIM_MODEL = 1: numerical solution, analytic theory and NCLASS \nSIM_MODEL = 2: numerical solution and analytic theory only \nSIM_MODEL = 3: analytic theory and NCLASS only. \nSIM_MODEL = 4: NN network of NEO DKE';

valeur.EQUILIBRIUM_MODEL	= 2;
type.EQUILIBRIUM_MODEL	    	= 'real';
borne.EQUILIBRIUM_MODEL		= {0,1,2,3};
defaut.EQUILIBRIUM_MODEL	= 2;
info.EQUILIBRIUM_MODEL		= 'Parameter which selects the geometric equilibrium model.\nEQUILIBRIUM_MODEL = 0: s-alpha \nEQUILIBRIUM_MODEL = 1: large aspect ratio \nEQUILIBRIUM_MODEL = 2: Miller \nEQUILIBRIUM_MODEL = 3: General Grad-Shafranov';

valeur.COLLISION_MODEL	= 4;
type.COLLISION_MODEL   	= 'real';
borne.COLLISION_MODEL	= {0,1,2,3,4,5};
defaut.COLLISION_MODEL	= 4;
info.COLLISION_MODEL	= 'Parameter which selects the collision operator model.\nCOLLISION_MODEL = 1: Connor model \nCOLLISION_MODEL = 2: zeroth-order Hirshman-Sigmar model \nCOLLISION_MODEL = 3: Full Hirshman-Sigmar model \nCOLLISION_MODEL = 4: Full linearized Fokker-Planck operator \nCOLLISION_MODEL = 5: FP test particle operator with ad hoc field particle operator';

valeur.N_ENERGY	= 6;
type.N_ENERGY  	= 'real';
borne.N_ENERGY  = [1,100];
defaut.N_ENERGY	= 6;
info.N_ENERGY	= 'The number of energy polynomials - 1 in the computational domain';

valeur.N_XI	= 17;
type.N_XI  	= 'real';
borne.N_XI      = [1,100];
defaut.N_XI	= 17;
info.N_XI	= 'The number of xi polynomials - 1 in the computational domain';

valeur.N_THETA	= 29;
type.N_THETA  	= 'real';
borne.N_THETA   = [1,100];
defaut.N_THETA	= 29;
info.N_THETA	= 'The number of theta gridpoints in the computational domain.';

valeur.DV	= 0;
type.DV	    	= 'real';
borne.DV		= {0,1,2};
defaut.DV	= 0;
info.DV		= 'Parameter which selects determines if and how D and V are computed. 0: only particle fluxes, 1: computes D and V assuming all impurities are traces, 2: linearly decompose flux/gradient relation for each impurity (correct method in general cases)';

interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet

param.valeur     = valeur;
param.type       = type;
param.borne      = borne;
param.defaut     = defaut;
param.info       = info;
param.interface  = interface;
param.description = 'GUI for NEO';   % description (une ligne) de la fonction
param.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
param.gui      ='';                             % nom de l'interface graphique specifique si elle existe
param.controle = '';                        % nom de la fonction de controle des valeurs si elle existe

