%Function used to define QLK-NN 10D parameters in METIS
function param = declare_QLKNN10D_parameters

% METIS implementation parameters
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

valeur.prediction_domain	= 'Restricted';
type.prediction_domain		= 'string';
borne.prediction_domain		= {'Full','Restricted','Extended'};
defaut.prediction_domain	= 'Restricted';
info.prediction_domain		= 'if = Full, use QLK-NN to predict from magnetic axis to x=0.9;\nif = Restricted , predict Te, Ti  and Ne only from x= 0.25 or sawtooth inversion radius if this one is more external to x= 0.9;\nif = Extended, predict from magnetic axis to x = 0.95 (last point before LCFS)';

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
   
valeur.coef_min		= 1e-3;
type.coef_min		= 'real';
borne.coef_min		= [0,10];
defaut.coef_min		= 1e-3;
info.coef_min		= 'minimum value of diffusion coefficient what ever is other tunning for convergence of the solver';
%label.coef_min

% Keep RAPTOR notation instead of previous one use in METIS to simplified user life
valeur.minChie	        = 0.1;
type.minChie		= 'real';
borne.minChie		= [0,1];
defaut.minChie		= 0.1;
info.minChie		= 'minimum value of Chi_e';

valeur.minChii	        = 0.1;
type.minChii		= 'real';
borne.minChii		= [0,1];
defaut.minChii		= 0.1;
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

% METIS do not used Di and Vi

% some information has been taken from RAPTOR interface to QLK-NN 10D
% use QLKNN-fortran default values for defaults. 
% and set parameters to values use by RAPTOR.
% some parameters must be recasted in logical before call of QLK-NN mexfile

valeur.apply_victor_rule	= 0;
type.apply_victor_rule		= 'real';
borne.apply_victor_rule		= {0,1};
defaut.apply_victor_rule	= 0;
info.apply_victor_rule		= 'Apply the ''Victor rule'' rotation model';
%label.apply_victor_rule

valeur.norms		= 1;
type.norms		= 'real';
borne.norms		= [0,1];
defaut.norms		= 1;
info.norms		= 'norms uses with ''Victor rule''';
%label.norms

valeur.calc_heat_transport	= 1;
type.calc_heat_transport	= 'real';
borne.calc_heat_transport	= {0,1};
defaut.calc_heat_transport	= 1;
info.calc_heat_transport	= 'Calculate q_e and q_i';
%label.calc_heat_transport

valeur.calc_part_transport	= 1;
type.calc_part_transport	= 'real';
borne.calc_part_transport	= {0,1};
defaut.calc_part_transport	= 1;
info.calc_part_transport	= 'Calculate Gamma/D/V';
%label.calc_part_transport

valeur.use_ETG	= 1;
type.use_ETG	= 'real';
borne.use_ETG	= {0,1};
defaut.use_ETG	= 1;
info.use_ETG	= 'Calculate ETG mode transport';
%label.use_ETG

valeur.use_ITG	= 1;
type.use_ITG	= 'real';
borne.use_ITG	= {0,1};
defaut.use_ITG	= 1;
info.use_ITG	= 'Calculate ITG mode transport';
%label.use_ITG

valeur.use_TEM	= 1;
type.use_TEM	= 'real';
borne.use_TEM	= {0,1};
defaut.use_TEM	= 1;
info.use_TEM	= 'Calculate TEM mode transport';
%label.use_TEM

valeur.constrain_inputs	= 1;
type.constrain_inputs	= 'real';
borne.constrain_inputs	= {0,1};
defaut.constrain_inputs	= 1;
info.constrain_inputs	= 'Clip inputs to min_input and max_input, bound by margin_input';
%label.constrain_inputs

valeur.constrain_outputs	= 1;
type.constrain_outputs		= 'real';
borne.constrain_outputs		= {0,1};
defaut.constrain_outputs	= 1;
info.constrain_outputs		= 'Clip outputs to min_input and max_input, bound by margin_input';
%label.constrain_outputs

valeur.margin_input	= 0.95;
type.margin_input	= 'real';
borne.margin_input	= [0,1];
defaut.margin_input	= 0.95;
info.margin_input	= 'Margin on input for clipping function';
%label.margin_input

valeur.min_output	= -100;
type.min_output		= 'real';
borne.min_output	= [-1000,1000];
defaut.min_output	= -100;
info.min_output		= 'Minimum value for output';
%label.min_output

valeur.max_output	= 100;
type.max_output		= 'real';
borne.max_output	= [-1000,1000];
defaut.max_output	= 100;
info.max_output		= 'Maximum value for output';
%label.max_output

valeur.margin_output	= 1;
type.margin_output	= 'real';
borne.margin_output	= [0,1];
defaut.margin_output	= 1;
info.margin_output	= 'Margin on ouput for clipping function';
%label.margin_output

valeur.em_stab		= 0;
type.em_stab		= 'real';
borne.em_stab		= {0,1};
defaut.em_stab		= 0;
info.em_stab		= 'Apply em stabilisation rule from JETTO. Currenty 0 (off) or 1 (reduce gradient by fast-ion pressure)';
%label.em_stab

valeur.merge_modes		= 1;
type.merge_modes		= 'real';
borne.merge_modes		= {0,1};
defaut.merge_modes		= 1;
info.merge_modes		= 'Merge ETG/ITG/TEM modes together to a single total flux (debugging parameter)';
%label.merge_modes

valeur.force_evaluate_all	= 0;
type.force_evaluate_all		= 'real';
borne.force_evaluate_all	= {0,1};
defaut.force_evaluate_all	= 0;
info.force_evaluate_all		= 'Force all networks to be evaluated (debugging parameter)';
%label.force_evaluate_all

valeur.min_input	= '[1.0, 0.0, 0.0, -5.0, 0.66, -1.0, 0.09, 0.25, -5.0, -100.0, -100.0]';
type.min_input		= 'string';
borne.min_input		= '';
defaut.min_input	= '[1.0, 0.0, 0.0, -5.0, 0.66, -1.0, 0.09, 0.25, -5.0, -100.0, -100.0]';
info.min_input		= 'Minimal values for each input';
%label.min_input

valeur.verbosity	= 0;
type.verbosityl		= 'real';
borne.verbosity		= {0,1,2,3,4,5,10};
defaut.verbosity	= 0;
info.verbosity		= 'Level of verbosity (debugging parameter)';
%label.verbosity

valeur.channel_tag	= 'QLK-NN 10D';
type.channel_tag	= 'string';
borne.channel_tag	= {'QLK-NN 10D'};
defaut.channel_tag	= 'QLK-NN 10D';
info.channel_tag	= 'Tag in plots for this network (debugging parameter)';
%label.channel_tag

valeur.max_input	= '[3.0, 14.0, 14.0, 6.0, 15.0, 5.0, 0.99, 2.5, 0.0, 100.0, 100.0]';
type.max_input		= 'string';
borne.max_input		= '';
defaut.max_input	= '[3.0, 14.0, 14.0, 6.0, 15.0, 5.0, 0.99, 2.5, 0.0, 100.0, 100.0]';
info.max_input		= 'Maximum values for each input';
%label.max_input

valeur.apply_stability_clipping		= 0;
type.apply_stability_clipping		= 'real';
borne.apply_stability_clipping		= {0,0};
defaut.apply_stability_clipping		= 0;
info.apply_stability_clipping		= 'currently dont care, pending removal';
%label.apply_stability_clipping

valeur.use_effective_diffusivity	= 0;
type.use_effective_diffusivity		= 'real';
borne.use_effective_diffusivity		= {0,0};
defaut.use_effective_diffusivity	= 0;
info.use_effective_diffusivity		= 'Use Gamma_e, not separate Ds and Vs';
%label.use_effective_diffusivity

valeur.factor_nustar	= 1;
type.factor_nustar	= 'real';
borne.factor_nustar	=  [0.1 10];
defaut.factor_nustar	= 1;
info.factor_nustar	= 'factor applied to nustar (to correct problem in collisionality operator for TEM): nustar -> nustar * factor_nustar';
%label.factor_nustar

valeur.factor_ETG	= 1;
type.factor_ETG	        = 'real';
borne.factor_ETG	=  [0 10];
defaut.factor_ETG	= 1;
info.factor_ETG	        = 'factor applied to ETG saturation rule (recommandation: 0.33)';
%label.factor_ETG

valeur.Effect_TioTe	    = 'real value';
type.Effect_TioTe	    = 'string';
borne.Effect_TioTe	    =  {'real value','forced to 2.5','forced to 0.25'};
defaut.Effect_TioTe	    = 'real value';
info.Effect_TioTe	    = 'This key is used to test effect of Ti/Te (in the limit of the NN training set) on transport coefficients;\nIn QLK-NN call, if not = real value, force the Ti/Te input to be equal to:\nif = forced to 2.5, Ti/Te = 2.5;if = forced to 0.25, Ti/Te = 0.25;\n';
%label.Effect_TioTe

valeur.use_ion_diffusivity_networks	= 0;
type.use_ion_diffusivity_networks	= 'real';
borne.use_ion_diffusivity_networks	= {0,0};
defaut.use_ion_diffusivity_networks	= 0;
info.use_ion_diffusivity_networks	= 'Use D_e and V_e, not D_i and D_e (if using Ds and Vs)';
%label.use_ion_diffusivity_networks

interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet

param.valeur     = valeur;
param.type       = type;
param.borne      = borne;
param.defaut     = defaut;
param.info       = info;
param.interface  = interface;
param.description = 'GUI for QLK-NN 10D';   % description (une ligne) de la fonction
param.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
param.gui      ='';                             % nom de l'interface graphique specifique si elle existe
param.controle = '';                        % nom de la fonction de controle des valeurs si elle existe

