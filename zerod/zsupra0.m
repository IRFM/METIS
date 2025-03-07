function [esupra,pth,tauseff] = zsupra0(temps,pinjo,taus,taue,ecrit,e0,agaz,flag_icrh,esupra_in)

  
if nargin < 8
   flag_icrh = 0;
end
if nargin < 9
   esupra_in = zeros(size(temps));
   evolution = 0;
else
   evolution = 1;
end

taus  = max(1e-3,taus);
taue  = max(1e-6,taue);
ecrit = max(30,ecrit);
e0    = max(1e3,e0);
agaz  = max(1,mean(agaz));

%  calcul de L. G . Ericksonn
mp          =   1.6726485e-27;            % masse au repos du proton (kg)
e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
vc          =   sqrt(2 .*e .*  ecrit ./ mp ./ agaz);
v0          =   sqrt(2 .*e .*  e0 ./ mp ./ agaz);
%tauseff_old     =   min(taus,100) .* (1 - (vc./v0) .^2 .* ...
%	       (1 ./ 3 .* log((v0 .^ 2 - v0 .* vc + vc .^2) ./ (v0 + vc) ./ vc) + ...
%		1 ./ sqrt(3) .* (atan( (2 .* v0 - vc) ./ vc ./ sqrt(3)) + atan(1./sqrt(3)))));
		
% nouvelle formule valide a toutes les energies.
% merci a Didier Moreau
x0 = max(eps,v0 ./ vc);
tauseff     =   min(taus,100) .* (1 + log((x0 .^ 2 + 2 .* x0  + 1) ./ (x0 .^ 2 - x0 + 1)) ./ 3 ./ x0 .^ 2 - ...
                2 .* (atan((2 .* x0 -1) ./ sqrt(3)) + atan(1./sqrt(3))) ./ sqrt(3) ./ x0 .^ 2);


%figure(21);clf; plot(x0,tauseff./min(taus,100),'.r',x0,tauseff_old./min(taus,100),'.b');drawnow;

% autre formule pour test (iter)
%tauseff     = taus ./ (1 + (vc./v0) .^ 3);

% regle empririque JET
if flag_icrh >0
   tauseff = tauseff ./ 4;
end
tauseff(tauseff <= 0)  = taue(tauseff <= 0) ./ 1e3;
tauref  = zpmean(temps,pinjo,tauseff);
tauseff = min(tauseff,10 .* tauref);
tauseff = (cat(1,tauref,tauseff(2:end))+tauseff) ./ 2;
tauseff(1) = tauref;
fon  = pinjo > (1e5 .* mean(taue));
tauseff(find(~fon)) = tauref;

% calcul uniquement si c'est interressant
if all(pinjo<eps)
   % pas de puissance
   pth = pinjo;
   esupra = zeros(size(pth));
   te = temps;
%  elseif all((tauseff.*1e2) < min(diff(temps)) )
%     pth = pinjo;
%     esupra = pth ./ 2 .* tauseff;
%     te = temps;
%     fprintf('S');
%     keyboard
else
   % integration pour le calcul de esupra
    esupra_ini  = esupra_in(1);
    if esupra_ini == 0
	  esupra_ini = 0.01 .* pinjo(1) ./2.*tauseff(1);
    end
   if  evolution == 1
	  esupra = esupra_in;
	  [te,esupra(2:end)] = z0ode(temps(2:end),pinjo(2:end),0.5 .* tauseff(2:end),esupra_in(2));
          esupra(1:2) = esupra_in(1:2);
	  emax        = cat(1,esupra_ini, esupra_in(2) + cumtrapz(temps(2:end),pinjo(2:end)));
   else
	  [te,esupra] = z0ode(temps,pinjo,0.5 .* tauseff,esupra_ini);
	  emax        = min(cumtrapz(temps,pinjo),max(pinjo) .* max(tauseff)/2);
	  esupra      = min(esupra,emax);
   end
   esupra      = min(esupra,emax);
   % puisssance
   % le 2 vient de l'integration
   %esupra    = pth .* tauseff ./ 2;
   % securite si non convergence
   if length(esupra) < length(pinjo)
	 fprintf('F');
         pth = pinjo;
         esupra = pinjo .* tauseff ./ 2;
   else
      pth          = esupra ./ tauseff .* 2;
   end
end

% securite
esupra(~isfinite(esupra)) = 0;
esupra(esupra<=0) = 0;
pth(~isfinite(pth))       = pinjo(~isfinite(pth));
pth                       = max(0,pth);
if  evolution == 0
	pth = min(max(pinjo),pth);
end

%  disp('in zupra0');
%  figure(21) ;clf;plot(temps,pth,temps,pinjo);drawnow
%  keyboard
%  function  s  = zpmean(t,p,e)
%  
%  indok = find(isfinite(p) & isfinite(e));
%  if isempty(indok)
%     s =1;
%  elseif length(indok) == 1
%  	s = e(indok);
%  else
%     s = trapz(t(indok),p(indok) .* e(indok)) ./ trapz(t(indok),eps + p(indok));
%  end
%  

function  s  = zpmean(t,p,e)

%  figure(21);clf
%  subplot(2,1,1)
%  plot(t,p)
%  subplot(2,1,2);
%  plot(t,e);
%  drawnow

indok = find(isfinite(p) & isfinite(e));
if length(indok) < 3
	indok = find(isfinite(e));
	if isempty(indok)
		s = NaN;
		disp('zpmean : invalid data ...');
	else
		s = mean(e(indok));
	end
else
	t = t(indok);
	p = p(indok);
	e = e(indok);
	indp = find(p > 0);
	%disp(length(indp))
	if length(indp) > 2
		s = trapz(t(indp),p(indp) .* e(indp)) ./ trapz(t(indp),eps + p(indp));
	else
		s = mean(e);
	end
end


%  function  s  = zpmean(t,p,e)
%  
%  indok = find(isfinite(p) & isfinite(e));
%  if isempty(indok)
%     s =1;
%  elseif length(indok) == 1
%  	s = e(indok);
%  else
%     s = trapz(t(indok),p(indok) .* e(indok)) ./ trapz(t(indok),eps + p(indok));
%  end
