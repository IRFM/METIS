% calcul un profil de q monotonic
function qp =z0qp(x,qmin,qa)

qmin = max(0.01,min(qmin,qa));

% Wesson p 114 (3.4)
if size(x,1) >1

   nup1 = qa ./ qmin;
   vt = ones(size(qa,1),1);
   ve = ones(1,size(x,2));
   qp = x .^ 2 ./ max(eps,(1 - (1 - x .^ 2) .^ (nup1 * ve)));
   qp =  ((qa ./max(qp(:,end),eps)) * ve) .* qp;
   qp(:,1) = qmin;

elseif min(x) < 0
   nup1 = qa ./ qmin;
   vt = ones(size(qa,1),1);
   ve = ones(1,size(x,2));
   qp = (vt * (x .^ 2)) ./ max(eps,(1 - (vt * (1 - x .^ 2)) .^ (nup1 * ve)));
   qp =  ((qa ./max(qp(:,end),eps)) * ve) .* qp;
   ind     = min(find(abs(x) == min(abs(x))));
   qp(:,ind) = qmin;

else
   nup1 = qa ./ qmin;
   vt = ones(size(qa,1),1);
   ve = ones(1,size(x,2));
   qp = (vt * (x .^ 2)) ./ max(eps,(1 - (vt * (1 - x .^ 2)) .^ (nup1 * ve)));
   qp =  ((qa ./max(qp(:,end),eps)) * ve) .* qp;
   qp(:,1) = qmin;
end  
