% Get NIST data via IAEA website : https://www-amdis.iaea.org/FLYCHK/
%FLYCHK calculations are known to give better results for highly ionized plasmas and intermediate electron densities.
%For very low Ne, FLYCHK may underestimate the radiative loss rates for near-neutral plasmas due to lack of ∆n=0 transitions.
%For very high Ne, FLYCHK may give inadequate results due to simple continuum lowering model and lack of degenerate states.
%Please note that FLYCHK uses a screened-hydrogenic level and results for near-neutral ions should be taken with caution. Especially for most elements with atomic numbers > 18, the ground state of neutral atom is incorrect and the near-neutral results are incorrect. 
%
%References :
% 1) FLYCHK: an extension to the K-shell spectroscopy kinetics model FLY, H. -K. Chung, W. L. Morgan and R. W. Lee, Journal of Quantitative Spectroscopy and Radiative Transfer, Volume 81, November 2003, Pages 107-115
% 2) FLYCHK: Generalized population kinetics and spectral model for rapid spectroscopic analysis for all elements, H.-K. Chung, M.H. Chen, W.L. Morgan, Y. Ralchenko and R.W. Lee, High Energy Density Physics, Volume 1, Issue 1, December 2005, Pages 3-12
%
function [Te,Z,rep] = get_iaea_nist_repartition(Z_element,density)

% init
Te     = [];
Z      = [];
rep    = [];

% density index
density_index = max(1,ceil(log10(density)) -6 - 11);

% URL
rep_url  = sprintf('https://www-amdis.iaea.org/FLYCHK/ZBAR/cp%3.3d.%2.2d.zvd',Z_element,density_index);
% get data
% temporary name
tpz = tempname;
% use wget to retreive the data
[s,t] = unix(sprintf('/usr/bin/wget -O %s %s; /usr/bin/cat %s; rm %s',tpz,rep_url,tpz,tpz));
if s~= 0
    error(t);
end

% loop on temperature
ind_temp = strfind(t,'fun: Te =');
indn  = strfind(t,newline);
indss = strfind(t,'//');
Te = NaN * ones(length(ind_temp),1);
rep = NaN * ones(Z_element+1,length(ind_temp));
for k=1:length(ind_temp)
    indl = ind_temp(k)+9;
    indf = min(indn(indn>indl));
    Te(k) = sscanf(t(indl:indf),'%g');
    %
    inddata1 = min(indss(indss > indf));
    l        = find(inddata1 == indss,1);
    inddata2 = strfind(t,sprintf('#end csd%d/u',k-1));
    %
    data1    = str2num(t((inddata1+3):(inddata2-2)));
    Z = data1(:,1);
    rep(:,k) = data1(:,2) ./ sum(data1(:,2));
end    
    


