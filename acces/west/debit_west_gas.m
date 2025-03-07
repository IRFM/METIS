function [time,flux_electrons,gastab] = debit_west_gas(Nch,temps)

% initialisation
gastab = [];
flux_electrons = [];
time = [];

% baraometry and gas command information
[bar,tbar]=tsbase(Nch,'GBARAT%8');
if ~isempty(bar)
    tbar=tbar(:,1);
else
    return
end

[vp,tvp]=tsbase(Nch,'GTENSION');
if isempty(vp)
  return
end
tvp=tvp(:,1);

[deb,tdeb]=tsbase(Nch,'GDEBIT');
if isempty(deb)
  return
end
tdeb=tdeb(:,1);
time = tdeb;
flux_electrons = zeros(size(time));

[dec,tdec]=tsbase(Nch,'GDEBCAL');
tdec=tdec(:,1);

[ptore,tptore]=tsbase(Nch,'GPTORE');
tptore=tptore(:,1);


% reading species
% line 1 to 11
l0111 = tsmat(Nch,'EXP=T=S;Fueling;Gas_D2');
l1216 = tsmat(Nch,'EXP=T=S;Fueling;Gas_IMP1');
indi  = find(l1216 > 0,1);
l1216 = l1216(indi);
l1721 = tsmat(Nch,'EXP=T=S;Fueling;Gas_IMP2');
indi  = find(l1721 > 0,1);
l1721 = l1721(indi);
lgas = cat(2,l0111 * ones(1,11),l1216 * ones(1,5),l1721 * ones(1,5));

%bar=bar-mean(bar(tbar>-10&tbar<-5));
%bar=smooth(bar,50);
ptore=ptore(:,1).*10.^(ptore(:,2)); ptore=smooth(ptore,20);

% physical constant
phys = cphys;

% loop on gas line
for i=1:size(vp,2)
    if (max(vp(:,i))>10.2) && (lgas(i) > 0)
        % Flow rate (Pa.m^3/s)
        [lg_out,Ze_out] = listgas(lgas(i));
        deb_loc = Ze_out .* deb(:,i) .* interp1(tvp,double(vp(:,i) > 10.2),time,'nearest',0);
        flux_electrons  = flux_electrons + deb_loc;
        if isfield(gastab,lg_out)
            gastab.(lg_out) = gastab.(lg_out) + phys.pam3 .* deb_loc;
        else
            gastab.(lg_out) = phys.pam3 .* deb_loc;
        end
    end
end
flux_electrons = phys.pam3 .* flux_electrons;

if nargin > 1
  indtime = find((time >= min(temps)) & time <= max(temps));
  if length(indtime) > 3
    width = max(3,ceil(mean(diff(temps)) /  mean(diff(time(indtime)))) .* 2 + 1);
  else
    width = 3;
  end
  %figure(21);clf;plot(time,flux_electrons,'b',time,sgolayfilt(flux_electrons,1,width),'r');drawnow
  flux_electrons = max(0,sgolayfilt(flux_electrons,1,width));
  flux_electrons = interp1(time,flux_electrons,temps,'nearest',0);
  noms = fieldnames(gastab);
  for k=1:length(noms)
      var = max(0,sgolayfilt(gastab.(noms{k}),1,width));
      gastab.(noms{k})= interp1(time,var,temps,'nearest',0);
  end
  time = temps;
end

function [lg_out,Ze_out] = listgas(index)

lg{1} = 'H2';
Ze(1) = 2;

lg{2} = 'D2';
Ze(2) = 2;

lg{3} = 'D2_0d1H2';  %'HE3';
Ze(3) = 2;

lg{4} = 'HE4';
Ze(4) = 2;

lg{5} = 'N2';
Ze(5) = 14;

lg{6} = 'Ne';
Ze(6) = 10;

lg{7} = 'Ar';
Ze(7) = 18;

lg{8} = 'Kr';
Ze(8) = 36;

lg{9} = 'Xe';
Ze(9) = 54;

lg{10} = 'CD4';
Ze(10) = 10;

if nargin == 0
    lg_out = lg;
    Ze_out = Ze;
else
    lg_out = lg{index};
    Ze_out = Ze(index);
end

%  strcpy(Nom_Gaz[ 0],"H2 \0");
%  strcpy(Nom_Gaz[ 1],"D2 \0");
%  strcpy(Nom_Gaz[ 2],"He3\0");
%  strcpy(Nom_Gaz[ 3],"He4\0");
%  strcpy(Nom_Gaz[ 4],"N2 \0");
%  strcpy(Nom_Gaz[ 5],"Ne \0");
%  strcpy(Nom_Gaz[ 6],"Ar \0");
%  strcpy(Nom_Gaz[ 7],"Kr \0");
%  strcpy(Nom_Gaz[ 8],"Xe \0");
%  strcpy(Nom_Gaz[ 9],"CD4\0")
%
%
function phys = cphys
% script d'estimation des flux de matiere
phys.c           =   2.99792458e8;             % speed of light in vacuum (m/s)  (definition)
phys.h           =   6.62606876e-34;           % Planck constant (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % electron charge (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeablity of vacuum (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivity of vacuum (F/m)  (definition)
phys.g           =   6.673e-11;                % gravitation constant (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % Boltzmann constant (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % fine structure constant (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % electron mass (at rest) (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % proton mass (at rest) (kg)
phys.ua          =   1.66053873e-27;           % Atomic mass unit (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % Avogadro number (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % Stephan constant ( W*m^-2*K^-4) (+/- 0.000040e-8)
phys.pam3        =   (4.41e-4 .* phys.avo);    % conversion d'un nombre de particules en en Pa.m^3

