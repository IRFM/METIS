function varliste = plotmhd1t(temps,varstr,postmhd,param)

% declaration des options
varliste = {'densite','vper','vpol','vpara','temp','Aper','Apol','Ator'};
if nargin < 4 
   return
end

% parametre de mishka
parametre = param.cons.mhd.stab;

% boucle sur les n
for k = 1:3
    if isfield(postmhd,sprintf('n%d',k));
        nk    = getfield(postmhd,sprintf('n%d',k));
        st =sprintf('choc %s simulation, a t = %g, nbrho = %d',param.from.machine, temps,size(nk.EV,3));
        set(gcf,'name',st);
	     l     = min(find(nk.t >= temps));
		  m     = strmatch(varstr,varliste,'exact');
		  parammish(1)  = parametre.nbrhostab;
		  parammish(4)  = parametre.rwall;
        parammish(5)  = parametre.ieta;
        parammish(6)  = parametre.ieq;
        parammish(7)  = parametre.dsurf1;
        parammish(8)  = parametre.dsurf;
        parammish(9)  = 0.8;
        parammish(10) = 0;
        parammish(11) = parametre.eta;
        parammish(12) = param.compo.z(1);
        parammish(2)  = parametre.polini(k);
        parammish(3)  = parametre.tor(k);

        ev  = squeeze(nk.EV(l,:,:));
        iev = squeeze(nk.IEV(l,:,:));
        subplot(2,3,k+3);
		  hold off
        plotharm1t(parammish,ev,iev,nk.EW(l),nk.IEW(l),m);
        subplot(2,3,k);
		  hold off
		  plot(nk.t,nk.EW,'r',nk.t,nk.IEW,'b')
		  hold on
		  xlabel('temps (s)');
		  ylabel('taux de croissance');
		  legend('Reel','Imag')
		  title(sprintf('N = %d, t =  %g, err = %g',parametre.tor(k),nk.t(l),nk.DEW(l)));
		  drawnow
	 else	 
	    subplot(2,3,k)
		 delete(gca); 
	    subplot(2,3,k+3)
		 delete(gca); 
	 end
end
