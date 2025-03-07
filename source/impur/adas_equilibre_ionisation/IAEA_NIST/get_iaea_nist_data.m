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
function [te,Zave,Lz] = get_iaea_nist_data(Z_element,density)

% init output empty 
te     = [];
Zave   = [];
Lz     = [];
% density is in cm^-3 instead of m^3
density = density ./ 1e6;

% URLs
Z_url  = sprintf('https://www-amdis.iaea.org/FLYCHK/ZBAR/zbar%3.3d.html',Z_element);
Z_data = read_data(Z_url);
% format data
indn = find(Z_data == max(Z_data));
Z_density = Z_data(1:indn);
nbc = indn +1;
nbl = length(Z_data((indn+1):end)) / nbc;
data_tab = reshape(Z_data((indn+1):end),nbc,nbl);
Z_te     = data_tab(1,:);
Zave_2d     = data_tab(2:end,:);
%
% interpolate on requested density
Zave = NaN * ones(1,length(Z_te));
for k=1:length(Z_te)
    Zave(k) = interp1(Z_density,Zave_2d(:,k),density,'pchip','extrap');
end

%
Lz_url = sprintf('https://www-amdis.iaea.org/FLYCHK/ZBAR/pwrt%3.3d.html',Z_element);
Lz_data = read_data(Lz_url);
% format data
indn = find(Z_data == max(Lz_data));
Lz_density = Lz_data(1:indn);
nbc = indn +1;
nbl = length(Lz_data((indn+1):end)) / nbc;
data_tab = reshape(Lz_data((indn+1):end),nbc,nbl);
Lz_te     = data_tab(1,:);
Lz_2d        = data_tab(2:end,:);
%
% interpolate on requested density
Lz = NaN * ones(1,length(Lz_te));
for k=1:length(Z_te)
    Lz(k) = interp1(Lz_density,Lz_2d(:,k),density,'pchip','extrap');
end
% convert from erg/atom/s to W*m^3
Lz =  1e-7  .* Lz ./ (density * 1e6);

if all(Z_te == Lz_te)
    te = Lz_te;
else
    error('data missmatch');
end



% read data
function data = read_data(turl)

% init
data =[];

% temporary name
tpz = tempname;
% use wget to retreive the data
[s,t] = unix(sprintf('/usr/bin/wget -O %s %s',tpz,turl));
if s~= 0
    error(t);
end
fid = fopen(tpz);
rep = fscanf(fid,'%s');
fclose(fid);
delete(tpz);
ind_1 = strfind(rep,'<tr><td>T<Sub>e</Sub>&#92;') + length('<tr><td>T<Sub>e</Sub>&#92;') + 1;
ind_end = strfind(rep,'</table>') -1;
rep = rep(ind_1:ind_end);

indf = strfind(rep,'<td>');
indp = strfind(rep,'</td>');
% loop on the field
data = [];
for k=1:length(indf)
    start = indf(k) + 4;
    stop = min(indp(indp > start))-1;
    if ~isempty(stop)
        loct = rep(start:stop);
        if ~isempty(loct)
            indsup = strfind(loct,'<Sup>');
            if isempty(indsup)
                %fprintf('%s =',loct);disp(sscanf(loct,'%g'))
                data(end+1) = sscanf(loct,'%g');
            else
                mant = str2double(loct(1:indsup-1));
                indnsup = strfind(loct,'</Sup>');
                expo  = str2double(loct((indsup+5):(indnsup-1)));
                %fprintf('%s =',loct);disp(mant .* 10 ^ expo)
                data(end+1) = mant .^ expo;
            end
        end
    end
end




