function choice = preset_jt60sa_scenario_choice

% init output
choice = 0;

if isdeployed && ~isempty(dir('reference_scenario_jt60-sa'))
    pathm = 'reference_scenario_jt60-sa';
elseif  ~isempty(dir('reference_scenario_jt60-sa'))
    pathm = 'reference_scenario_jt60-sa';
elseif ~isempty(dir(fullfile(fileparts(which('jt60sawall.mat')),'reference_scenario_jt60-sa')))
    pathm = fullfile(fileparts(which('jt60sawall.mat')),'reference_scenario_jt60-sa');
else
    % nodata
    return
end

list     = {'1','2','3','4-1','4-2','5-1','5-2','6','last setting'};
list_num = {1,2,3,4.1,4.2,5.1,5.2,6};
choice   = menu('Which reference scenario you want simulate ?',list{:});
switch choice 
case {0,length(list)}
    % cancel case
    return
end

filename = fullfile(pathm,sprintf('METIS_JT-60SA_scenario_%s.mat',list{choice}));
if ~exist(filename,'file')
    error(sprintf('unable to find %s',filename));
end
info = load(filename);
if ~isfield(info,'post')
  error(sprintf('%s is not a METIS file !',filename));
end
post = info.post;


% already set parameters
if evalin('base','exist(''jt60sa_param'')')
  jt60sa_param = evalin('base','jt60sa_param');
else
  jt60sa_param =[];
end

% parameters data set source for METIS
valeur.reference_parameters = filename;
% deafault set of parameters for selected scenario
valeur.sepa_option = list_num{choice};      
if ~isempty(jt60sa_param) && (list_num{choice} == 2)
    valeur.sepa_create = jt60sa_param.sepa_create;    
else
    valeur.sepa_create = 'No';    
end
valeur.gas = post.z0dinput.option.gaz;   % gas species as in METIS (1=H, 2=D, 3=DT & 4=He)      
valeur.ip     = round(max(post.z0dinput.cons.ip) / 1e4) / 1e2;   % plasma current (MA)  
switch valeur.sepa_option
% choix de la forme du plasma en plateau
case 1
  b0 = 2.25;
case 2
  b0 = 2.25;
case 3
  b0 = 2.25;
case 4.1
  b0 = 2.28;
case 4.2
  b0 = 2.28;
case 5.1
  b0 = 1.72;             
case 5.2
  b0 = 1.62;
case 6
  b0 = 1.41;
otherwise
  error('unknown scenario');
end
valeur.b0     = b0; 
if ~isempty(jt60sa_param)
    valeur.voltage_first = jt60sa_param.voltage_first;    
else
    valeur.voltage_first = 'high';  
end
if valeur.ip >= 2.5
    valeur.flux = 'full';    
else
    valeur.flux = 'half';    
end
valeur.density     = sum(post.z0dinput.cons.nbar .* (post.z0dinput.cons.ip == max(post.z0dinput.cons.ip))) ./ ...
		     sum((post.z0dinput.cons.ip == max(post.z0dinput.cons.ip)));  
valeur.density     = round(valeur.density ./ 1e17) / 100;		     
valeur.edge_density_factor     = round(post.z0dinput.option.nea_factor * 100) ./ 100;   
valeur.H_H     = round(max(post.z0dinput.cons.hmore .* (post.z0dinput.cons.ip == max(post.z0dinput.cons.ip))) *100) ./ 100; 
if post.z0dinput.option.sitb ~= 0
    valeur.ITB = 'on';   
else
    valeur.ITB = 'off';   
end
valeur.PNBI_N     = round(max(real(post.z0dinput.cons.pnbi) .* (post.z0dinput.cons.ip == max(post.z0dinput.cons.ip))) / 1e5) / 10;   
valeur.PNBI_P     = round(max(imag(post.z0dinput.cons.pnbi) .* (post.z0dinput.cons.ip == max(post.z0dinput.cons.ip))) / 1e5) / 10;   
valeur.PICRH      = round(max(post.z0dinput.cons.picrh .* (post.z0dinput.cons.ip == max(post.z0dinput.cons.ip))) / 1e5) / 10;   
valeur.PECCD      = round(max(post.z0dinput.cons.pecrh .* (post.z0dinput.cons.ip == max(post.z0dinput.cons.ip))) / 1e5) / 10;   
valeur.PBREAK     = round(post.z0dinput.cons.pecrh(1) / 1e5) /10;
dd                = abs(post.z0dinput.cons.ip - 0.8 * max(post.z0dinput.cons.ip));
indramp           = find(dd == min(dd),1) - 1;
valeur.PRAMPUP    = round(post.z0dinput.cons.pecrh(indramp) / 1e5) / 10;   
valeur.Recycling  = round(post.z0dinput.option.Recycling * 1000) / 1000; 
if post.z0dinput.option.matthews  ==  1
  valeur.radiation = 'Matthews';
else
  valeur.radiation = 'Lz';
end
valeur.SOL_model = post.z0dinput.option.sol_model; 
%finally always of by default
if post.z0dinput.option.runaway == 0
      valeur.runaway = 'off';  
else
      valeur.runaway = 'on'; 
end
%valeur.runaway = 'off';  
%finally always of by default
if post.z0dinput.option.berror == 0
    valeur.breakdown = 'off';   
else
    valeur.breakdown = 'on';   
end
%valeur.breakdown = 'off';   
valeur.duration     = round(post.z0dinput.cons.temps(max(find((post.z0dinput.cons.ip == max(post.z0dinput.cons.ip))))+1));
valeur.f_dipdt_rampup     = 1;   
valeur.f_dipdt_rampdown     = 1;   

zassignin('base','jt60sa_param',valeur);
