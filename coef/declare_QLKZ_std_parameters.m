%Function used to define QLK-NN 10D parameters in METIS
function param = declare_QLKZ_std_parameters

% METIS implementation parameters

valeur.QLKZ_mode	= 'Full';
type.QLKZ_mode		= 'string';
borne.QLKZ_mode		= {'Full','NN-like'};
defaut.QLKZ_mode	= 'Full';
info.QLKZ_mode		= 'if = Full, use QLKZ for all METIS species; \n = NN-like, use QLKZ with 10D neural network setup with only one impurity (carbon) whose content is set from Zeff';

valeur.smooth_space	= 1;
type.smooth_space	= 'real';
borne.smooth_space	= {0,1};
defaut.smooth_space	= 1;
info.smooth_space	= 'if = 1, applied a filter on space direction to smooth transport coeffecient';

valeur.anti_aliasing		= 0.1;
type.anti_aliasing		= 'real';
borne.anti_aliasing		= [0,1];
defaut.anti_aliasing		= 0.1;
info.anti_aliasing		= 'anti aliasing factor';
%label.anti_aliasing
  
valeur.temperature_factor	= 1e-6;
type.temperature_factor		= 'real';
borne.temperature_factor	= [0,1];
defaut.temperature_factor	= 1e-6;
info.temperature_factor		= 'temperature factor for simulated annealing';
%label.temperature_factor
  
valeur.nbloop_max	= 501;
type.nbloop_max		= 'real';
borne.nbloop_max	= [0,1];
defaut.nbloop_max	= 501;
info.nbloop_max		= 'maximum number of loop';
%label.nbloop_max

valeur.prediction_domain	= 'Full';
type.prediction_domain		= 'string';
%borne.prediction_domain		= {'Full','Restricted','Clamped'};
borne.prediction_domain		= {'Full','Restricted'};
defaut.prediction_domain	= 'Full';
info.prediction_domain		= 'if = Full, use QLKZ to predict up to magnetic axis otherwise in restricted mode, predict Te, Ti  and Ne only up to x= 0.25 or sawtooth inversion radius if this one is mode external';

valeur.domain_rho		= 1;
type.domain_rho	    	        = 'real';
borne.domain_rho		= [0,1];
defaut.domain_rho		= 1;
info.domain_rho		        = 'Normalised radial coordinate for the boundary condition. Must be used with prediction domain=Full';

valeur.factor_neo_chii		= 1;
type.factor_neo_chii		= 'real';
borne.factor_neo_chii		= [0,100];
defaut.factor_neo_chii		= 1;
info.factor_neo_chii		= 'factor applied ot Chi_neo (computed with Hinton formula) that is added to turbulent coefficient Chi_ion';
%label.factor_neo_chii

valeur.neo_ion_in_D		= 0;
type.neo_ion_in_D		= 'real';
borne.neo_ion_in_D		= [0,1];
defaut.neo_ion_in_D		= 0;
info.neo_ion_in_D		= 'fraction of Chi_ion_neo limiting the minimal value of De';
%label.neo_ion_in_D

valeur.neo_ion_in_Chie		= 0;
type.neo_ion_in_Chie		= 'real';
borne.neo_ion_in_Chie		= [0,1];
defaut.neo_ion_in_Chie		= 0;
info.neo_ion_in_Chie		= 'fraction of Chi_ion_neo limiting the minimal value Chi_e';
%label.neo_ion_in_Chie
   
valeur.minChie	        = 0.05;
type.minChie		= 'real';
borne.minChie		= [0,1];
defaut.minChie		= 0.05;
info.minChie		= 'minimum value of Chi_e';

valeur.minChii	        = 0.05;
type.minChii		= 'real';
borne.minChii		= [0,1];
defaut.minChii		= 0.05;
info.minChii		= 'minimum value of Chi_i';

valeur.minDe	        = 0.05;
type.minDe		= 'real';
borne.minDe		= [0,1];
defaut.minDe		= 0.05;
info.minDe		= 'minimum value of D_e';

valeur.minVe	        = -100;
type.minVe		= 'real';
borne.minVe		= [-1000,0];
defaut.minVe		= -100;
info.minVe		= 'minimum value of V_e';

valeur.maxChie	        = 100;
type.maxChie		= 'real';
borne.maxChie		= [1,1000];
defaut.maxChie		= 100;
info.maxChie		= 'maximum value of Chi_e';

valeur.maxChii	        = 100;
type.maxChii		= 'real';
borne.maxChii		= [1,1000];
defaut.maxChii		= 100;
info.maxChii		= 'maximum value of Chi_i';

valeur.maxDe	        = 100;
type.maxDe		= 'real';
borne.maxDe		= [1,1000];
defaut.maxDe		= 100;
info.maxDe		= 'maximum value of D_e';

valeur.maxVe	        = 100;
type.maxVe		= 'real';
borne.maxVe		= [1,1000];
defaut.maxVe		= 100;
info.maxVe		= 'maximum value of V_e';


valeur.coef_min		= 1e-3;
type.coef_min		= 'real';
borne.coef_min		= [0,10];
defaut.coef_min		= 1e-3;
info.coef_min		= 'minimum value of diffusion coefficient what ever is other tunning for convergence of the solver';
%label.coef_min

valeur.calc_heat_transport	= 1;
type.calc_heat_transport	= 'real';
borne.calc_heat_transport	= {0,1};
defaut.calc_heat_transport	= 1;
info.calc_heat_transport	= 'Calculate q_e and q_i';
%label.calc_heat_transport

valeur.calc_part_transport	= 0;
type.calc_part_transport	= 'real';
borne.calc_part_transport	= {0,1};
defaut.calc_part_transport	= 0;
info.calc_part_transport	= 'Calculate Gamma/D/V';
%label.calc_part_transport

valeur.ETG_mult	= 1;
type.ETG_mult	= 'real';
borne.ETG_mult	= [-100,100];
defaut.ETG_mult	= 1;
info.ETG_mult	= 'ETG flux multiplier';
%label.ETG_mult

valeur.em_stab		= 0;
type.em_stab		= 'real';
borne.em_stab		= {0,1};
defaut.em_stab		= 0;
info.em_stab		= 'Apply em stabilisation rule from JETTO. Currenty 0 (off) or 1 (reduce gradient by fast-ion pressure)';
%label.em_stab

valeur.verbosity	= 0;
type.verbosityl		= 'real';
borne.verbosity		= {0,1,2,3,4,5,10};
defaut.verbosity	= 0;
info.verbosity		= 'Level of verbosity (debugging parameter)';
%label.verbosity

valeur.use_effective_diffusivity	= 0;
type.use_effective_diffusivity		= 'real';
borne.use_effective_diffusivity		= {0,0};
defaut.use_effective_diffusivity	= 0;
info.use_effective_diffusivity		= 'Use Gamma_e, not separate Ds and Vs';
%label.use_effective_diffusivity

valeur.coll_mult	= 1;
type.coll_mult  	= 'real';
borne.coll_mult	        =  [0,100];
defaut.coll_mult	= 1;
info.coll_mult	        = 'factor applied to collisionality';
%label.coll_mult

valeur.channel_tag	= 'QLKZ';
type.channel_tag	= 'string';
borne.channel_tag	= {'QLKZ'};
defaut.channel_tag	= 'QLKZ';
info.channel_tag	= 'Tag in plots for this model';
%label.channel_tag


valeur.Effect_TioTe	= 'real value';
type.Effect_TioTe	    = 'string';
borne.Effect_TioTe	    =  {'real value','forced to 0.25','forced to 1'};
defaut.Effect_TioTe	    = 'real value';
info.Effect_TioTe	    = 'This key is used to test effect of Ti/Te (in the limit of the NN training set) on transport coefficients;\nIn QLK-NN call, if not = real value, force the Ti/Te input to be equal to:\nif = forced to 0.25, Ti/Te = 0.25;if = forced to 1, Ti/Te = 1;\n';
%label.Effect_TioTe

valeur.domain_TiovTe	    = 0.7;
type.domain_TiovTe	    = 'real';
borne.domain_TiovTe	    =  [0,1];
defaut.domain_TiovTe	    = 0.7;
info.domain_TiovTe	    = 'Choose the max for which Ti/Te is modified in the transport model';
%label.domain_TiovTe

valeur.relacc1		= 1e-3;
type.relacc1	    	= 'real';
borne.relacc1		= [0,1];
defaut.relacc1		= 1e-3;
info.relacc1		= 'Relative accuracy in 1D integrals';

valeur.relacc2		= 2e-2;
type.relacc2	    	= 'real';
borne.relacc2		= [0,1];
defaut.relacc2		= 2e-2;
info.relacc2		= 'Relative accuracy in 2D integrals. Major factor in setting runtime!';

valeur.numsols		= 3;
type.numsols	    	= 'real';
borne.numsols		= {1,2,3};
defaut.numsols		= 3;
info.numsols		= 'Maximum number of unstable solutions per wavenumber in output. \nTypically no more than 2 instabilities appear at a given wavenumber, since we are limited to electrostatic instabilities in QuaLiKiz (e.g. ITG and TEM)';

valeur.phys_meth	= 2;
type.phys_meth	    	= 'real';
borne.phys_meth		= {0,1,2};
defaut.phys_meth	= 2;
info.phys_meth		= 'If 0, outputs only total fluxes. \nIf 1, outputs additional separate particle diffusivity and pinches. \nIf 2, further outputs separate heat diffusivity and pinches.';

valeur.coll_flag	= 1;
type.coll_flag	    	= 'real';
borne.coll_flag		= {0,1};
defaut.coll_flag	= 1;
info.coll_flag		= 'If 0, collisionless simulation (slight speed-up due to 1D electron trapped integral). \nIf 1, collisional simulation';

valeur.rot_flag	        = 0;
type.rot_flag	    	= 'real';
borne.rot_flag		= {0,1,2};
defaut.rot_flag	        = 0;
info.rot_flag		= 'WARNING, only available option=0 \nIf 0, rotation not included in dispersion relation (but still included for heavy impurity asymmetries). \nIf 1, rotation is fully included. \nIf 2, then rotation is only included for the outer half radius. \nThis is default for integrated modelling due to the known limitation of QuaLiKiz to resolve the parallel-velocity-gradient destabilisation term, more important at inner radii (see J. Citrin et al., PPCF 2017).';

valeur.separateflux	= 1;
type.separateflux    	= 'real';
borne.separateflux	= {0,1};
defaut.separateflux	= 1;
info.separateflux	= 'Output ITG/TEM seperated fluxes. See output. ETG separated flux is also provided by default';

valeur.numproc  	= 8;
type.numproc    	= 'real';
borne.numproc   	= [0,100];
defaut.numproc  	= 8;
info.numproc    	= 'Number of processors to be used for QuaLiKiz simulations';



interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet

param.valeur     = valeur;
param.type       = type;
param.borne      = borne;
param.defaut     = defaut;
param.info       = info;
param.interface  = interface;
param.description = 'GUI for QuaLiKiz';   % description (une ligne) de la fonction
param.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
param.gui      ='';                             % nom de l'interface graphique specifique si elle existe
param.controle = '';                        % nom de la fonction de controle des valeurs si elle existe

