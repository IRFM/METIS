% fonction auxiliaire pour le calul de nhem
function dndt = ztf0(t,n,void,tv,sv,tauv)

ind    = min(find(t<=tv));
if isempty(ind)
   ind = length(tv);
end
s      = sv(ind);
tau    = tauv(ind);
dndt   = - n ./ tau + s;

%if t <= 1.5
%  keyboard
%end
