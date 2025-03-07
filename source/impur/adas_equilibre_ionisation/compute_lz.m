% computation of ionisation equilibrium
% input:
%  element = element name (ex: 'Li')
%  year    = year of the computing methode (leave empty for default)
%  source  = source of the data: 'adas' for ADAS data base, 'web' for open ADAS and 'mattioli' for Mattioli data (default = mattioli, if empty auto select database). 
%  nH0ne   = neutral H densty over electron density (default = 0)
%  density = plasma density (default = 1e19, small effect) (m^-3)
%  Te      = vector of temperature points for sampling (if empty used = logspace(-3,5,201)) (eV)
%
% output: 
%  Te      = vector of temperature points for sampling (if empty used = logspace(-3,5,201))
%  Lz      = radaited power coefficient or radiative cooling rate
%  Zave    = averaged charge
%  Z2ave   = averaged squared charge
%  Zl      = state chage from 0 to number of charge of selected element
%  rep     = matrix of fractionnal abundances for each temperature and state charge
%  year    = methode used for primary data computation 
%  source  = selected source for data
%  quality = quality measurement
function [Te,Lz,Zave,Z2ave,Zl,rep,year_out,source_out,quality] = compute_lz(element,year,source,nH0ne,density,Te)

% sortie
Lz    = [];
Zave  = [];
Z2ave = [];


% gestion des entrees
if nargin < 2
  year = [];
end
if nargin < 3
  source = 'mattioli';
end
if nargin < 4
  nH0ne = 0;
end
if nargin < 5
  density = 1e19 / 1e6;
else
  density = density / 1e6;
end
% vecteur de temperature
if (nargin < 6) || isempty(Te)
  Te = logspace(-3,5,201);
end

% coputation of equilibrium
if length(year) > 1
  if year(1) == 0;
      year_loc = [];
  else
      year_loc = year(1);
  end
  if year(2) == 0;
      year = [];
  else
      year = year(2);
  end
else
  year_loc = year;
end
[Te,Zl,rep,year_out,source_out,quality] = equi_ionisation(element,year_loc,source,nH0ne,density,Te);
if isempty(rep)
       source_out = {source_out};
       year_out   = {year_out};
       return
end
% reading collisional radiative coefficients
[Te_plt,Dens,plt,HeaderComment,year_out2,source_out2] = read_adf11(element,'plt',year,source);
% extraction de la densite (petit effet)
ind_dens = find(Dens >= density,1);
if isempty(ind_dens)
    ind_dens = length(Dens);
end
plt = squeeze(plt(ind_dens,:,:));
if any(size(plt) == 1);
    plt = plt';
end
[Te_prb,Dens,prb,HeaderComment,year_out3,source_out3] = read_adf11(element,'prb',year,source);
ind_dens = find(Dens >= density,1);
if isempty(ind_dens)
    ind_dens = length(Dens);
end
prb = squeeze(prb(ind_dens,:,:));
if any(size(prb) == 1);
    prb = prb';
end
[Te_prc,Dens,prc,HeaderComment,year_out4,source_out4] = read_adf11(element,'prc',year,source);
if isempty(prc)
  prc = zeros(size(prb));
  Te_prc = Te_prb;
else
    ind_dens = find(Dens >= density,1);
    if isempty(ind_dens)
	ind_dens = length(Dens);
    end
    prc = squeeze(prc(ind_dens,:,:));
    if any(size(prc) == 1);
	prc = prc';
    end
end
source_out = {source_out,source_out,source_out2,source_out3,source_out4};
year_out   = {year_out,year_out,year_out2,year_out3,year_out4};

% preparation interpolation
plt    = cat(1,max(0,interp1(log10(Te_plt),plt,-16,'linear','extrap')),plt);   
Te_plt = union(eps,Te_plt(:));
if max(Te_plt) < 1e5
      plt = cat(1,plt,max(0,interp1(log10(Te_plt),plt,5,'linear','extrap')));   
      Te_plt = union(Te_plt(:),1e5);
end
prb    = cat(1,max(0,interp1(log10(Te_prb),prb,-16,'linear','extrap')),prb);   
Te_prb = union(eps,Te_prb(:));
if max(Te_prb) < 1e5
      prb = cat(1,prb,max(0,interp1(log10(Te_prb),prb,5,'linear','extrap')));   
      Te_prb = union(Te_prb(:),1e5);
end
prc    = cat(1,max(0,interp1(log10(Te_prc),prc,-16,'linear','extrap')),prc);   
Te_prc= union(eps,Te_prc(:));
if max(Te_prc) < 1e5
      prc = cat(1,prc,max(0,interp1(log10(Te_prc),prc,5,'linear','extrap')));   
      Te_prc = union(Te_prc(:),1e5);
end

% interpolation
Te = Te(:);
plt    = interp1(log10(Te_plt),plt,log10(Te),'pchip','extrap');
prb    = interp1(log10(Te_prb),prb,log10(Te),'pchip','extrap');
prc    = interp1(log10(Te_prc),prc,log10(Te),'pchip','extrap');

% combinaison
rad_coef = zeros(length(Te),size(plt,2) + 1);
rad_coef(:,1:(end-1))  =  plt;
rad_coef(:,2:end)      =  rad_coef(:,2:end) + prb;
if numel(nH0ne) == 1
    rad_coef(:,2:end)      =  rad_coef(:,2:end) + nH0ne .* prc;
else
     rad_coef(:,2:end)     =  rad_coef(:,2:end) + (nH0ne * ones(1,size(prc,2))) .* prc;   
end

%figure(22);clf;loglog(Te_plt,plt,'r',Te_prb,prb,'b',Te_prc,prc,'k');drawnow
%figure(22);clf;semilogx(Te_plt,plt,'r',Te_prb,prb,'b',Te_prc,prc,'k');drawnow

% compute output
vt      = ones(length(Te),1);
zmat    = vt * Zl;
Lz      = sum(rad_coef .* rep,2);
Zave    = sum(zmat .* rep,2);
Z2ave   = sum(zmat .^ 2 .* rep,2);


