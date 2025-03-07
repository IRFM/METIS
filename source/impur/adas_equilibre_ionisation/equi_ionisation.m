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
%  Zl      = state chage from 0 to number of charge of selected element
%  rep     = matrix of fractionnal abundances for each temperature and state charge
%  year    = methode used for primary data computation 
%  source  = selected source for data
%  quality = quality measurement


function [Te,Zl,rep,year_out,source,quality] = equi_ionisation(element,year,source,nH0ne,density,Te)

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
quality = [];

switch source
case 'mattioli'
 % charge de l'element
 [A,Z,name] = chargemasse(element);
 try
    if Z < 18
	rep = equicoronal_leger(Z,Te,1);     
    else
	rep = equicoronal_lourd(Z,Te,1);
    end
    ind_ok = find(all(isfinite(rep),2));
    rep    = rep(ind_ok,:);
    Te     = Te(ind_ok);
    quality = NaN .* ones(size(Te));
    year_out = 0;
 catch	
      rep = [];
      quality = NaN .* ones(size(Te));
      year_out = NaN;        
      Zl = 0:(size(rep,2) - 1);
 end
otherwise
    % lecture des tables
    [Te_i,Dens,ionisation,HeaderComment,year_out,source] = read_adf11(element,'scd',year,source);
    if isempty(ionisation)
        quality = NaN .* ones(size(Te));
        year_out = NaN;
        rep = [];
        Zl = 0:(size(rep,2) - 1);
        return
    end
    % extraction de la densite (petit effet)
    ind_dens = find(Dens >= density,1);
    if isempty(ind_dens)
	ind_dens = length(Dens);
    end
    ionisation = squeeze(ionisation(ind_dens,:,:));
    if any(size(ionisation) == 1);
	ionisation = ionisation';
    end
    [Te_r,Dens,recombinaison,HeaderComment,year_out,source] = read_adf11(element,'acd',year,source);
    % extraction de la densite (petit effet)
    ind_dens = find(Dens >= density,1);
    if isempty(ind_dens)
	ind_dens = length(Dens);
    end
    recombinaison = squeeze(recombinaison(ind_dens,:,:));
    if any(size(recombinaison) == 1);
	recombinaison = recombinaison';
    end
    [Te_x,Dens,xechange,HeaderComment,year_out_void,source] = read_adf11(element,'ccd',year,source);
    if isempty(xechange)
	xechange = zeros(size(recombinaison));
	Te_x = Te_r;
    else
	% extraction de la densite (petit effet)
	ind_dens = find(Dens >= density,1);
	if isempty(ind_dens)
	    ind_dens = length(Dens);
	end
	xechange = squeeze(xechange(ind_dens,:,:));
	if any(size(xechange) == 1);
	    xechange= xechange';
	end
    end
    v0       = zeros(1,size(ionisation,2));
    v1       = ones(1,size(ionisation,2));
    %Te   = union(eps,union(union(Te_i(:),Te_r(:)),Te_x(:)));
    %Te   = logspace(-3,5,201);
    if max(Te) < 1e5
      Te = union(Te,1e5);
    end
    %ionisation = cat(1,ionisation(1,:),ionisation);
    ionisation = cat(1,v0,ionisation);
    Te_i = union(eps,Te_i(:));
    if max(Te_i) < 1e5
      ionisation = cat(1,ionisation,max(0,interp1(log10(Te_i),ionisation,5,'linear','extrap')));   
      Te_i = union(Te_i(:),1e5);
      %ionisation = cat(1,ionisation,ionisation(end,:));   
      %ionisation = cat(1,ionisation,v0);   
    end
    %recombinaison = cat(1,recombinaison(1,:),recombinaison);   
    recombinaison = cat(1,max(0,interp1(log10(Te_r),recombinaison,-16,'linear','extrap')),recombinaison);   
    Te_r = union(eps,Te_r(:));
    if max(Te_r) < 1e5
      Te_r = union(Te_r(:),1e5);
      %recombinaison = cat(1,recombinaison,recombinaison(end,:));   
      recombinaison = cat(1,recombinaison,v0);   
    end
    Te_x = union(eps,Te_x(:));
    %xechange = cat(1,xechange(1,:),xechange);
    xechange = cat(1,v0,xechange);
    if max(Te_x) < 1e5
      xechange = cat(1,xechange,max(0,interp1(log10(Te_x),xechange,5,'linear','extrap')));   
      Te_x = union(Te_x(:),1e5);
    end
    %figure(21);semilogx(Te_i,ionisation,'b',Te_r,recombinaison,'r',Te_x,xechange,'k');drawnow
    % interpolation
    Te = Te(:);
    ionisation    = interp1(log10(Te_i),ionisation,log10(Te),'pchip','extrap');
    recombinaison = interp1(log10(Te_r),recombinaison,log10(Te),'pchip','extrap');
    xechange      = interp1(log10(Te_x),xechange,log10(Te),'pchip','extrap');
    % melange recombinaison et charge echange
    if numel(nH0ne) == 1
	  xrec = recombinaison + nH0ne .* xechange;
    else
	  xrec = recombinaison + (nH0ne * ones(1,size(xechange,2))) .* xechange;
	  
    end
    % boucle sur les temperature
    nc = size(ionisation,2) + 1;
    rep = NaN * ones(length(Te),nc);
    for k=1:length(Te)
      % matrice 
      mat_equi = zeros(nc,nc);
      for l=1:nc
	  if l <= 1
	      mat_equi(l,l)   = - ionisation(k,l);
	      mat_equi(l,l+1) =   xrec(k,l); 
	  elseif l >= nc
	      mat_equi(l,l-1) =    ionisation(k,l-1);
	      mat_equi(l,l)   =  - xrec(k,l-1);
	  else
	      mat_equi(l,l-1) =   ionisation(k,l-1);
	      mat_equi(l,l)   = - ionisation(k,l) - xrec(k,l-1);
	      mat_equi(l,l+1) =   xrec(k,l); 
	  end
      end
      %sol = null(mat_equi ./ max(abs(mat_equi(:))));
      %if all(size(sol)>1)
      %	  keyboard
      %end  
      [V,e]= eig(mat_equi./ max(abs(mat_equi(:))));
      ind_ok = find(diag(e) == 0);
      if isempty(ind_ok)
	sol = null(mat_equi./ max(abs(mat_equi(:))));
	qs  = sum(sol,1);
	indbest = max(find(abs(qs -1) == min(abs(qs -1))));
	sol = sol(:,indbest);
	if sol(max(abs(sol)) == abs(sol)) < 0
	  sol = - sol;
	end
      else
	sol = V(:,ind_ok);
	qs  = sum(sol,1);
	indbest = find(abs(qs -1) == min(abs(qs -1)),1);
	sol = sol(:,indbest);
	if sol(max(abs(sol)) == abs(sol)) < 0
	  sol = - sol;
	end      
      end
      % final normalisation
      quality(k) = abs(sum(sol) - 1);
      sol = sol' ./ sum(sol);
      sol((sol < eps) & (sol ~= 0)) = eps;
      sol = sol ./ sum(sol);
      rep(k,:) = sol;
    end
end
Zl = 0:(size(rep,2) - 1);


