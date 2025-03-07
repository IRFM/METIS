function choice = preset_west_scenario_choice

% init output
choice = 0;

if isdeployed && ~isempty(dir('reference_scenario_WEST'))
    pathm = 'reference_scenario_WEST';
elseif  ~isempty(dir('reference_scenario_WEST'))
    pathm = 'reference_scenario_WEST';
elseif ~isempty(dir(fullfile(fileparts(which('west_wall.mat')),'reference_scenario_WEST')))
    pathm = fullfile(fileparts(which('west_wall.mat')),'reference_scenario_WEST');
else
    % nodata
    return
end

ld = dir(fullfile(pathm,'WEST_METIS_*.mat'));
list =[];
for k =1:length(ld)
  list{k} = ld(k).name; 
end
list{end+1} = 'last setting';
list_num = 1:length(list);
rep = dir(fullfile(fileparts(which('westmetissimulation.m')),'configuration','*.txt'));
conf_list=[];
for kl =1:length(rep)
  conf_list{kl}  = rep(kl).name;
end
conf_list{end+1} = 'default tuning of the scenario generator';
choice   = menu('Which reference scenario you want simulate ?',list{:});
switch choice 
case {0,length(list)}
    % cancel case
    return
end

filename = fullfile(pathm,list{choice});
if ~exist(filename,'file')
    error(sprintf('unable to find %s',filename));
end
info = load(filename);
if ~isfield(info,'post')
  error(sprintf('%s is not a METIS file !',filename));
end
post = info.post;


% already set parameters
if evalin('base','exist(''west_param'')')
  west_param = evalin('base','west_param');
else
  info = westmetissimulation;
  west_param = info.valeur;
end


% debut creation structure
valeur.LCFS = west_param.LCFS;
valeur.configuration = conf_list{choice};
valeur.gas = post.z0dinput.option.gaz;   % gas species as in METIS (1=H, 2=D, 3=DT & 4=He)      
valeur.ip     = round(max(post.z0dinput.cons.ip) / 1e4) / 1e2;   % plasma current (MA)  
valeur.b0     = round(trapz(post.z0dinput.cons.temps,post.z0dinput.cons.ip .* post.z0dinput.geo.b0 .* post.z0dinput.geo.R) ./  ...
                trapz(post.z0dinput.cons.temps,post.z0dinput.cons.ip  .* post.z0dinput.geo.R) .* 10) ./ 10; 
valeur.ip_first = west_param.ip_first;
valeur.f_Greenwald = west_param.f_Greenwald;
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
valeur.PICRH      = round(max(post.z0dinput.cons.picrh .* (post.z0dinput.cons.ip == max(post.z0dinput.cons.ip))) / 1e5) / 10;   
valeur.PLHCD      = round(max(post.z0dinput.cons.plh .* (post.z0dinput.cons.ip == max(post.z0dinput.cons.ip))) / 1e5) / 10;   
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
%finally always of by default
if post.z0dinput.option.berror == 0
    valeur.breakdown = 'off';   
else
    valeur.breakdown = 'on';   
end
valeur.duration     = round(post.z0dinput.cons.temps(max(find((post.z0dinput.cons.ip == max(post.z0dinput.cons.ip))))+1));
valeur.Zeff         = max(post.z0dinput.cons.zeff .* (post.z0dinput.cons.ip == max(post.z0dinput.cons.ip))); 
valeur.Cw           = post.z0dinput.option.cw_offset;

valeur.f_dipdt_rampup    = west_param.f_dipdt_rampup;   
valeur.f_dipdt_rampdown  = west_param.f_dipdt_rampdown;   

zassignin('base','west_param',valeur);


