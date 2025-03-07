function [qnew,indmix] = z0dds(x,psi,q,mode,w1,epsq)

% si pas de dds
psinew = psi;
qnew   = q;
indmix = 0;

% calcul du flux helicoidal :
psistar=cumtrapz(psi,q - 1);

% Calcul de la position de q=1
[psistarmax,ind1]=max(psistar);
%figure(21);clf;plot(psi,psistar);drawnow
ind1    = max(ind1);
x1      = x(ind1);
aa = find(psistar<0);
if isempty(aa)
	% pas de dds
	return
end
aa = aa(aa>ind1);   % to avoid taking too small mixing lentgh indmix, in case of non monotonic q-profiles
if isempty(aa)
	% pas de dds
	return
end
indmix = aa(1);  % F.I. 13/12/2004
%indmix  = min(find(psistar<0)); % old version 
xmix    = x(indmix);
% Valeurs par defaut
indout = indmix - 1;
indin =1;

if ind1==1
	% pas de dds
	return
elseif x1 < 0.1
   % pas de DDS (trop petit)
   return
end

% Porcelli model
if mode == 1
   d         = abs(psistar -psistar(1));
   d(1:ind1) = inf;

   % calcul de la position de la zone de shear nul
   psicons    = (1-w1) .* psistarmax;
   dg         = abs(psistar(1:ind1) -psicons);
   indin      = max(find(min(dg) == dg));
   xin        = x(indin);
   indm       = ind1:indmix;
   dd         = abs(psistar(indm) -psicons);
   indout     = min(indm(min(find(min(dd) == dd))),indmix-1);
   xout       = x(indout);
   % profil de q
   if indin > 1
     nq0  = xin .^ 4 ./ 4 ./ trapz(x(1:indin),x(1:indin) .^ 3 ,2);
     q0   = xin .^ 4 ./ 4 ./ trapz(x(1:indin),x(1:indin) .^ 3 ./ q(1:indin),2)./nq0;
   else
     q0   = 1. ;
   end  
   dq        = epsq .* linspace(-0.5,0.5,length(indin:indout)); 
   
   if indin == 1
          xq        = cat(2,x(indin:indout),x(indmix:end));
          qq        = cat(2,ones(size(indin:indout)) + dq,q(indmix:end));
          qnew      = pchip(xq,qq,x);
   elseif indin == 2
          xq        = cat(2,x(1),x(indin:indout),x(indmix:end));
          qq        = cat(2,q0,ones(size(indin:indout)) + dq,q(indmix:end));
          qnew      = pchip(xq,qq,x);
   else                    
          xq        = cat(2,x(1:(indin-1)),x(indin:indout),x(indmix:end));
          qq        = cat(2,q0 .*ones(1,indin-1),ones(size(indin:indout)) + dq,q(indmix:end));
          qnew      = pchip(xq,qq,x);
   end

elseif mode == 2
   % Kadomtsev model
   xxint=linspace(0,x1,50);
   psistarint=interp1(x,psistar,xxint,'linear','extrap');
   xxext=interp1(psistar(ind1:end),x(ind1:end),psistarint,'linear','extrap');
   xxf=sqrt(xxext.^2-xxint.^2);

   psistarf = psistar;
   % il faudrait imposer gradient(psistarf,x)=0 au centre
   psistarf(1:indmix-1) = interp1(xxf,psistarint,x(1:indmix-1),'linear','extrap');

   qf     = 1. + pdederive(psi,psistarf,0,2,2,1);
   qnew   = q;
   qnew(1:indmix-1) = qf(1:indmix -1);
   
else
    % Incomplete Kadomtsev model
   % calcul de la position de la zone Kadomtsev
   psicons    = (1-w1) .* psistarmax;
   dg         = abs(psistar(1:ind1) -psicons);
   indin      = max(find(min(dg) == dg));
   xin        = x(indin);
   indm       = ind1:indmix;
   dd         = abs(psistar(indm) -psicons);
   indout     = min(indm(min(find(min(dd) == dd))),indmix-1);
   xout       = x(indout);
   % Zone Taylor : q = q0
   if indin > 1
     nq0       = xin .^ 4 ./ 4 ./ trapz(x(1:indin),x(1:indin) .^ 3 ,2);
     q0        = xin .^ 4 ./ 4 ./ trapz(x(1:indin),x(1:indin) .^ 3 ./ q(1:indin),2)./nq0;
   else
     q0        = 1. ;
   end
   
   % Zone Kadomtsev
   xxint=linspace(xin,x1,100);
   psistarint=interp1(x,psistar,xxint,'linear','extrap');
   xxext=interp1(psistar(ind1:end),x(ind1:end),psistarint,'linear','extrap');
   xxf=sqrt(xin^2+xxext.^2-xxint.^2);
   
   psistarf=psistar;
   psistarf(indin+1:indout-1) = interp1(xxf,psistarint,x(indin+1:indout-1),'linear','extrap');
   psistarf(indin)            = psistarmax ;
   
   if indin > 2
     psistartaylor = cumtrapz(psi(1:indin),q0*ones(size(x(1:indin)))-1.);
     psistarf(1:indin-1) = cumtrapz(psi(1:indin-1),q0*ones(size(x(1:indin-1)))-1.);
     % Continuity between Taylor relaxed zone and Kadomtsev zone
     psistarf(1:indin-1) = psistartaylor(1:indin-1)+psistarf(indin)-psistartaylor(indin);
   end
   
   qf     = 1. + pdederive(psi,psistarf,0,2,2,1) ;
   qnew   = q ;
   if indin > 1
     qnew(1:indin-1) = q0*ones(size(x(1:indin-1))) ;
   end
   qnew(indin)            = 1. ;
   qnew(indin+1:indmix-1) = qf(indin+1:indmix -1);
end

%  inte   = - F .* vpr .* r2i ./ (4*pi^2) ./ qnew;
%  psinew = rhomax .* cumtrapz(x,inte,2);
%  % raccordement a rmix :
%  psinew(1:indmix)       = psinew(1:indmix) - psinew(indmix) + psi(indmix);
%  psinew((indmix+1):end) = psi((indmix+1):end);

