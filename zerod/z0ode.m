% integartion rapide de dn/dt = - n/tau + s en mode discret
% test de la fonction :
% t= linspace(0,100,10001)';
% s= cos(t) .* exp(sin(t));
% tau = 0.1 + pi .* (s >0);
% [t0,r0] = z0ode(t,s,tau,pi,0);
% [t1,r1] = z0ode(t,s,tau,pi,1);
% [t2,r2] = z0ode(t,s,tau,pi,2);
% figure;clf;plot(t0,r0,'r',t1,r1,'b',t2,r2,'k')
function [trep,rep] = z0ode(temps,source,tau,ini,speed)

if nargin < 5
   speed = 1;
end

if speed == -1
   dt   = diff(temps);
   ek1   = exp(- dt ./ tau(1:(end-1)));
   ek2   = exp(- dt ./ tau(2:end));
   rep  = ini(1) .* ones(size(temps));
   for  k = 1:(length(temps)-1)
         % shema explicite
         rep1 = rep(k) .* ek1(k) + source(k) .* tau(k) .* (1 - ek1(k));
         % shema implicite
         rep2 = rep(k) .* ek2(k) + source(k+1) .* tau(k+1) .* (1 - ek2(k));
         % melange
         rep(k+1) = (rep1 + rep2) ./ 2;
    end 
    trep = temps;

elseif speed == 1
   dt   = diff(temps);
   ek1   = min(1,exp(- dt ./ tau(1:(end-1))));
   ek2   = min(1,exp(- dt ./ tau(2:end)));
   rep  = ini(1) .* ones(size(temps));
   for  k = 1:(length(temps)-1)
         % shema explicite
         rep1 = rep(k) .* ek1(k) + source(k) .* tau(k) .* (1 - ek1(k));
         % shema implicite
         rep2 = rep(k) .* ek2(k) + source(k+1) .* tau(k+1) .* (1 - ek2(k));
         % melange
         rep(k+1) = (rep1 + rep2) ./ 2;
    end 
    trep = temps;
elseif speed == 2
   dt   = diff(temps);
   rep  = ini(1) .* ones(size(temps));
   n1   = ini(1);
   tau1 = tau(1);
   s1   = source(1);
   for  k = 1:(length(temps)-1)
         tau0 = tau1;
         tau1 = tau(k+1);
         s0   = s1;
         s1   = source(k+1);
         n0   = n1;
         n1 = -(tau0+tau1.*dt(k)).^(-1/tau1).*tau0^(1/tau1).*(s0.*tau0+2.*s0.*tau0.*tau1-s1.*tau0.^2-2.*n0.*tau1.^2-3.*n0.*tau1-n0)./(2.*tau1.^2+3.*tau1+1)+ ...
               (2.*dt(k).*s0.*tau1.^2+s1.*dt(k).^2.*tau1.^2+s0.*tau1.*dt(k)+s1.*dt(k).^2.*tau1+2.*s0.*tau0.*tau1+s0.*tau0-s1.*tau0.^2+s1.*dt(k).*tau0)./(2.*tau1+1)./(tau1+1);
         rep(k+1) = n1;        
   end 
   trep = temps;

else
   [trep,rep] = ode23('ztf0',temps,ini(1),[],temps,source,tau);
   if length(rep) < length(temps)
      dt   = diff(temps);
      ek1   = min(1,exp(- dt ./ tau(1:(end-1))));
      ek2   = min(1,exp(- dt ./ tau(2:end)));
      rep  = ini(1) .* ones(size(temps));
      for  k = 1:(length(temps)-1);
         % shema explicite
         rep1 = rep(k) .* ek1(k) + source(k) .* tau(k) .* (1 - ek1(k));
         % shema implicite
         rep2 = rep(k) .* ek2(k) + source(k+1) .* tau(k+1) .* (1 - ek2(k));
         % melange
         rep(k+1) = (rep1 + rep2) ./ 2;
      end 
    trep = temps;
   end
end


% pour test
return
figure(19)
clf
[tc,rc] = ode23('ztf0',temps,ini(1),[],temps,source,tau);
plot(trep,rep,'r',tc,rc,'b',temps,source.*tau,'g');
drawnow
if length(rc) == length(rep)
   if (abs(mean(rc - rep))./ abs(mean(rc))) > 0.1
      keyboard
   end
end
