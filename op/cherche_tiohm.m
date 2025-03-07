function [ne,te,ti0,tim,ti1] = cherche_tiohm(fr)

vloop = 1;
r0    = 2.34;
epar  = vloop /(2*pi*r0)

% boucle sur ne
for  k = 1:10
	nec = k *1e19;
	
	% boucle sur te
	for l=1:10
		tec = 500 *l;
		tic = linspace(500,tec);
		% pour le centre
		[qei,qohm,eta] = comp_source(0,nec,tec,tic,epar);
		indok          = find(qei <(fr*qohm));
		if ~isempty(indok)
			ti0(k,l) = min(tic(indok));
		else
			ti0(k,l) =NaN;
		end
		% pour le milieu
		[qei,qohm,eta] = comp_source(0.5,nec,tec,tic,epar);
		indok          = find(qei <(fr*qohm));
		if ~isempty(indok)
			tim(k,l) = min(tic(indok));
		else
			tim(k,l) =NaN;
		end
		% pour le bord
		[qei,qohm,eta] = comp_source(1,nec,tec,tic,epar);
		indok          = find(qei <(fr*qohm));
		if ~isempty(indok)
			ti1(k,l) = min(tic(indok));
		else
			ti1(k,l) =NaN;
		end
		
		ne(k,l) = nec;
		te(k,l) = tec;
	end
end
