function ti = zborneti(x,te,ne,ti,epar)

% par defaut
fr =1;


% malheureusement 2 boucles
for k =1:size(ti,1)
	ec =max(abs(epar(k)),mean(abs(epar)));
	for l =1:size(ti,2)
		nec = max(ne(k,l),1e18);
		tec = max(te(k,l),100);
		xc  = x(1,l);
		tic = max(min(te(k,:)),13.6):100:tec;
	   [qei,qohm,eta] = comp_source(xc,nec,tec,tic,ec);
	   indok          = find(qei <(fr*qohm));
		if ~isempty(indok)
			ti(k,l) =  max(13.6,max(min(tec,ti(k,l)),min(tic(indok))));
		else
			ti(k,l) =  max(13.6,min(tec,ti(k,l)));
		end
	end
	ti(k,:)=zmonotone(x(1,:),ti(k,:));
end
