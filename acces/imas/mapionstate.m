function state = mapionstate(Z,A,label,te,ti,nions,nfast,psupra,psupra_perp,psupra_para,vtor,vpol,state)

if nargin < 13
    warning('Sub structure ion.state is not initialized');
    state = {};
end

% constante physique (phys)
phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)

% initialisation
persistent tabmat
% nouvelle donnees ADAS
if isempty(tabmat)
   if exist('Lz_zave.mat','file')
       load('Lz_zave.mat','tabmat')
   else
       load('abundances.mat','tabmat')   
   end
end


% search if element is in tabmat
noms = fieldnames(tabmat);
ind  = min(strmatch(label,noms,'exact'));
stateone = state{1};
if ~isempty(state)
    velocity = state{1}.velocity;
end
if isempty(ind)
  fprintf('unable to find element %s in data\n',label);
  % put NaN data
  stateone.z_min =  1;
  stateone.z_max =  Z;
  stateone.z_average= -9.0000e+40;
  stateone.z_square_average =  -9.0000e+40;
  stateone.ionisation_potential = -9.0000e+40;
  stateone.label = label;
  stateone.name = label;
  stateone.electron_configuration = [];
  stateone.vibrational_level = -9.0000e+40;
  stateone.vibrational_mode = [];
  velocity.radial= [];
  velocity.rdiamagnetic = [];
  velocity.parallel     = [];
  velocity.poloidal     = vpol;
  velocity.toroidal     = vtor; 
  stateone.velocity     = velocity;
  stateone.temperature  = ti;
  stateone.density      = nions + nfast;
  stateone.density_thermal = nions;
  stateone.density_fast    = nfast;
  stateone.pressure        = phys.e .* nions .* ti + psupra;
  stateone.pressure_thermal = phys.e .* nions .* ti;
  stateone.pressure_fast_perpendicular = psupra_perp;
  stateone.pressure_fast_parallel = psupra_para; 
  state{1} = stateone;

else
  % loop on charge states without neutral
  info = tabmat.(noms{ind});
  tek  = te ./ 1e3;
  frac = NaN * ones(length(tek),length(info.Zl));
  for k=1:length(info.Zl)
    if info.Zl(k) == 0
      frac(:,k) = 0; % only plasma 
    else
      frac(:,k) = interp1(log10(info.Te),info.rep(:,k),log10(tek),'pchip','extrap');
    end
  end
  % normalisation
  nf = sum(frac,2);
  %disp([mean(nf),std(nf)]);
  frac = frac ./ (nf * ones(1,size(frac,2)));
  % filling
  for k=2:length(info.Zl)
      if ~isempty(state)
        stateone = state{1};
        velocity = state{1}.velocity;
      end
      stateone.z_min =  info.Zl(k);
      stateone.z_max =  info.Zl(k);
      stateone.z_average= info.Zl(k);
      stateone.z_square_average =  info.Zl(k) .^ 2;
      stateone.ionisation_potential = -9.0000e+40;
      stateone.label = label;
      stateone.name = label;
      stateone.electron_configuration = [];
      stateone.vibrational_level = -9.0000e+40;
      stateone.vibrational_mode = [];
      velocity.radial= [];
      velocity.rdiamagnetic = [];
      velocity.parallel     = [];
      velocity.poloidal     = vpol .* frac(:,k)';
      velocity.toroidal     = vtor .* frac(:,k)'; 
      stateone.velocity     = velocity;
      stateone.temperature  = ti;
      stateone.density      = (nions + nfast) .* frac(:,k)';
      stateone.density_thermal = nions .* frac(:,k)';
      stateone.density_fast    = nfast .* frac(:,k)';
      stateone.pressure        = phys.e .* nions .* ti.* frac(:,k)' + psupra .* frac(:,k)';
      stateone.pressure_thermal = phys.e .* nions .* ti .* frac(:,k)';
      stateone.pressure_fast_perpendicular = psupra_perp .* frac(:,k)';
      stateone.pressure_fast_parallel = psupra_para .* frac(:,k)';
      state{end+1} = stateone;
  end
    
  
  
end

%size(state)
%state{1}.density'