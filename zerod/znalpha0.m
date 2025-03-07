% calcul le contenu en helium du plamsa
function na = znalpha0(t,s,tau,nini,vp,evolution)

if nargin < 6
    evolution = 0;
end

n0      =  nini(1) .* vp(1);
% limite  sur tau
tau   = max(1e-3,abs(tau(:)));
%tau     = max(0.1,tau);
s     = max(1,s);
%s       = max(max(1,max(s)./1e5),s);
%[ta,na_alt] = ode23t('ztf0',t,n0,[],t,s,tau);
[ta,na] = z0ode(t,s,tau,n0);

%na - na_alt
%disp('in znalpha0')
%keyboard
