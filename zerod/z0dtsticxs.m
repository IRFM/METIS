% script de plot de la temperature ionique
shot  = fix(post.z0dinput.shot);
[tic,stic] = cgcgettrait(fix(shot),'ticxs');
if isempty(tic.times)
	
	ticxs = load(sprintf('/home/sccp/gcbd/fenzi/CXRS/ks4fit_java/analyses/TICXS_%d',fix(shot)));
	if any(~isfinite(ticxs.pos))
		ticxs.pos = (ones(size(ticxs.Ticxs,1),1) * [2.36 2.45 2.53 2.61 2.67 2.74 2.80 2.87]);
		ticxs.erpos = 0.2 .* ones(size(ticxs.pos));
	end
elseif isfield(tic,'rhofit')

	ticxs.pos       = tic.rhofit;
	ticxs.erpos     = 0.02 .* ones(size(ticxs.pos)); 
	ticxs.Ticxs     = tic.ticxs;
	ticxs.erTicxs   = tic.erticxs;
	ticxs.tcxs      = tic.times;
	ticxs.fluxC     = tic.fluxcxs;
	ticxs.erfluxC   = tic.erfluxcxs;
	ticxs.wphi      = -tic.wphicxs;
	ticxs.erwphi    = tic.erwphicxs;
	if isempty(stic.raycxs.temps)
		ticxs.pos = (ones(size(ticxs.Ticxs,1),1) * [2.45 2.53 2.61 2.67 2.74 2.80 2.87]);
		ticxs.erpos = 0.2 .* ones(size(ticxs.pos));	
	end
	
else
	ticxs.pos       = ones(size(tic.times)) * stic.raycxs.temps';
	ticxs.erpos     = ones(size(tic.times)) * stic.erraycxs.temps'; 
	ticxs.Ticxs     = tic.ticxs;
	ticxs.erTicxs   = tic.erticxs;
	ticxs.tcxs      = tic.times;
	ticxs.fluxC     = tic.fluxcxs;
	ticxs.erfluxC   = tic.erfluxcxs;
	ticxs.wphi      = -tic.wphicxs;
	ticxs.erwphi    = tic.erwphicxs;
	if isempty(stic.raycxs.temps)
		ticxs.pos = (ones(size(ticxs.Ticxs,1),1) * [2.45 2.53 2.61 2.67 2.74 2.80 2.87]);
		ticxs.erpos = 0.2 .* ones(size(ticxs.pos));	
	end
	
end

if ~isempty(ticxs)
	
	if size(ticxs.pos,1)==1
	    ticxs.pos       = ones(size(ticxs.tcxs)) *  ticxs.pos;
	    ticxs.erpos     = ones(size(ticxs.tcxs)) *  ticxs.erpos; 	
	end
	
	% recherche des indices dans la simulation cronos
	indl = (1:length(post.profil0d.temps))';
	indc = interp1(post.profil0d.temps,indl,ticxs.tcxs,'nearest');
	indl = (1:length(post.zerod.temps))';
	indz = interp1(post.zerod.temps,indl,ticxs.tcxs,'nearest');
	% creation de la figure (une par temps)	
	k = 1;
	h = findobj(0,'type','figure','tag',sprintf('zticxsti%d',k));
	if isempty(h)
		h=figure('tag',sprintf('zticxsti%d',k));
	else
		figure(h);
	end
	clf
	set(h,'defaultaxesfontsize',16,'defaultaxesfontweight', ...
	      'bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
		
	m = max(1,round(sqrt(length(indc))));
	n = max(1,ceil(length(indc) ./ m));
	ligne = 1;
	col   = 1;
	for k = 1:length(indc)
		subplot(m,n,k)
		% changement de coordonnees     
		amat = interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'nearest');   
		Rti = post.profil0d.Raxe(indc(k),:) + amat(indc(k)) .* post.profil0d.xli;
		plot(Rti,post.profil0d.tip(indc(k),:),'b',ticxs.pos(k,:),ticxs.Ticxs(k,:),'or');
		if k == 1
		  legend('Metis','Ticxs');
		end
		hold on
		Rti = post.profil0d.Raxe(indc(k),:) - amat(indc(k)) .* post.profil0d.xli;
		plot(Rti,post.profil0d.tip(indc(k),:),'b');
		% creation de vecteur pour le plot des erreur
        	tion=[ticxs.Ticxs(k,:);ticxs.Ticxs(k,:);ones(size(ticxs.Ticxs(k,:)))*NaN];
        	tion=tion(:);
        	pos=[ticxs.pos(k,:)+ticxs.erpos(k,:);ticxs.pos(k,:)-ticxs.erpos(k,:);ones(size(ticxs.pos(k,:)))*NaN];
        	pos=pos(:);
		plot(pos,tion,'r');
        	tion=[ticxs.Ticxs(k,:)+ticxs.erTicxs(k,:);ticxs.Ticxs(k,:)-ticxs.erTicxs(k,:);ones(size(ticxs.Ticxs(k,:)))*NaN];
        	tion=tion(:);
        	pos=[ticxs.pos(k,:);ticxs.pos(k,:);ones(size(ticxs.pos(k,:)))*NaN];
        	pos=pos(:);
		plot(pos,tion,'r');
		
		if ligne == m
		    xlabel('R (m)');
		end
		if col == 1
		    ylabel('T_i (eV)');
		end
		if k == 1
		    title(sprintf('Tore Supra, shot # %d @ %g s',fix(shot),post.profil0d.temps(indc(k))));
		else
		    title(sprintf('@ %g s',post.profil0d.temps(indc(k))));
		end
		set(gca,'ButtonDownFcn','zdataplot(''extrait'');');
		
		col = col + 1;
		if col > n
		  col = 1;
		  ligne = ligne +1;
		end
		
        end
	drawnow
	
	% plot de la rotation
	k = 1;
	h = findobj(0,'type','figure','tag',sprintf('zticxsrot%d',k));
	if isempty(h)
		h=figure('tag',sprintf('zticxsrot%d',k));
	else
		figure(h);
	end
	clf
	set(h,'defaultaxesfontsize',16,'defaultaxesfontweight', ...
	      'bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
	% creation de la figure (une par temps)	
	ligne = 1;
	col   = 1;
	for k = 1:length(indc)
		subplot(m,n,k)
		% cahngement de coordonnees     
		amat = interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'nearest');   
		Rti = post.profil0d.Raxe(indc(k),:) + amat(indc(k)) .* post.profil0d.xli;
		plot(Rti,post.profil0d.vtor(indc(k),:)./Rti ,'b', ...
		ticxs.pos(k,:),ticxs.wphi(k,:),'or',Rti,post.profil0d.omega(indc(k),:),'g');
		if k == 1
		  legend('Metis','Ticxs');
		end
		hold on
		Rti = post.profil0d.Raxe(indc(k),:) - amat(indc(k)) .* post.profil0d.xli;
	        plot(Rti,post.profil0d.vtor(indc(k),:)./Rti,'.b',Rti,post.profil0d.omega(indc(k),:),'.g');
		% creation de vecteur pour le plot des erreur
        	tion=[ticxs.wphi(k,:);ticxs.wphi(k,:);ones(size(ticxs.wphi(k,:)))*NaN];
        	tion=tion(:);
        	pos=[ticxs.pos(k,:)+ticxs.erpos(k,:);ticxs.pos(k,:)-ticxs.erpos(k,:);ones(size(ticxs.pos(k,:)))*NaN];
        	pos=pos(:);
		plot(pos,tion,'r');
        	tion=[ticxs.wphi(k,:)+ticxs.erwphi(k,:);ticxs.wphi(k,:)-ticxs.erwphi(k,:);ones(size(ticxs.wphi(k,:)))*NaN];
        	tion=tion(:);
        	pos=[ticxs.pos(k,:);ticxs.pos(k,:);ones(size(ticxs.pos(k,:)))*NaN];
        	pos=pos(:);
		plot(pos,tion,'r');
		
		plot(post.profil0d.Raxe(indc(k),1),post.zerod.wrad(indz(k)),'*c');
		
		if ligne == m
		    xlabel('R (m)');
		end
		if col == 1
		    ylabel('Wphi (m/s)');
		end
		if k == 1
		  title(sprintf('Tore Supra, shot # %d @ %g s',fix(shot),post.profil0d.temps(indc(k))));
		else
		    title(sprintf('@ %g s',post.profil0d.temps(indc(k))));
		end
		set(gca,'ButtonDownFcn','zdataplot(''extrait'');');
		
		col = col + 1;
		if col > n
		  col = 1;
		  ligne = ligne +1;
		end
		      
        end
        drawnow
        
	% plot de la rotation
%  	k = 2;
%  	h = findobj(0,'type','figure','tag',sprintf('zticxsrot%d',k));
%  	if isempty(h)
%  		h=figure('tag',sprintf('zticxsrot%d',k));
%  	else
%  		figure(h);
%  	end
%  	clf
%  	set(h,'defaultaxesfontsize',16,'defaultaxesfontweight', ...
%  	      'bold','defaultaxesfontname','times', ...
%  	      'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
%  	% creation de la figure (une par temps)	
%  	ligne = 1;
%  	col   = 1;
%  	for k = 1:length(indc)
%  		subplot(m,n,k)
%  		% cahngement de coordonnees     
%  		amat = interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'nearest');   
%  		Rti = post.profil0d.Raxe(indc(k),:) + amat(indc(k)) .* post.profil0d.xli;
%  		plot(Rti,post.profil0d.vtor(indc(k),:)./Rti ,'b', ...
%  		ticxs.pos(k,:),-ticxs.wphi(k,:),'or',Rti,post.profil0d.omega(indc(k),:),'g');
%  		if k == 1
%  		  legend('Metis','Ticxs');
%  		end
%  		hold on
%  		Rti = post.profil0d.Raxe(indc(k),:) - amat(indc(k)) .* post.profil0d.xli;
%  	        plot(Rti,post.profil0d.vtor(indc(k),:)./Rti,'.b',Rti,post.profil0d.omega(indc(k),:),'.g');
%  		% creation de vecteur pour le plot des erreur
%          	tion=[ticxs.wphi(k,:);ticxs.wphi(k,:);ones(size(ticxs.wphi(k,:)))*NaN];
%          	tion=tion(:);
%          	pos=[ticxs.pos(k,:)+ticxs.erpos(k,:);ticxs.pos(k,:)-ticxs.erpos(k,:);ones(size(ticxs.pos(k,:)))*NaN];
%          	pos=pos(:);
%  		plot(pos,-tion,'r');
%          	tion=[ticxs.wphi(k,:)+ticxs.erwphi(k,:);ticxs.wphi(k,:)-ticxs.erwphi(k,:);ones(size(ticxs.wphi(k,:)))*NaN];
%          	tion=tion(:);
%          	pos=[ticxs.pos(k,:);ticxs.pos(k,:);ones(size(ticxs.pos(k,:)))*NaN];
%          	pos=pos(:);
%  		plot(pos,-tion,'r');
%  		
%  		plot(post.profil0d.Raxe(indc(k),1),post.zerod.wrad(indz(k)),'*c');
%  		
%  		if ligne == m
%  		    xlabel('R (m)');
%  		end
%  		if col == 1
%  		    ylabel('Wphi (m/s)');
%  		end
%  		if k == 1
%  		    title(sprintf('Tore Supra, shot # %d @ %g s',fix(shot),post.profil0d.temps(indc(k))));
%  		else
%  		    title(sprintf('@ %g s',post.profil0d.temps(indc(k))));
%  		end
%  		set(gca,'ButtonDownFcn','zdataplot(''extrait'');');
%  		
%  		col = col + 1;
%  		if col > n
%  		  col = 1;
%  		  ligne = ligne +1;
%  		end
%  		      
%          end
%          drawnow

	% plot pour le carbonne
	k=1;
	h = findobj(0,'type','figure','tag',sprintf('zticxsc%d',k));
	if isempty(h)
		h=figure('tag',sprintf('zticxsc%d',k));
	else
		figure(h);
	end
	clf
	set(h,'defaultaxesfontsize',16,'defaultaxesfontweight', ...
	      'bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
	ligne = 1;
	col   = 1;
		
	% creation de la figure (une par temps)	
	for k = 1:length(indc)
		subplot(m,n,k)
		% cahngement de coordonnees     
		amat = interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'nearest');   
		Rti = post.profil0d.Raxe(indc(k),:) + amat(indc(k)) .* post.profil0d.xli;
		if post.z0dinput.option.zimp == 6
		    fluxc = post.profil0d.nzp(indc(k),:)./ (6 .* post.profil0d.nep(indc(k),:)) .*  ...
			    abs(post.profil0d.s0(indc(k),:)+post.profil0d.s0m(indc(k),:));
		else
		    fluxc = post.z0dinput.option.rimp .* post.profil0d.nzp(indc(k),:)./ (6 .* post.profil0d.nep(indc(k),:)) .*  ...
			    abs(post.profil0d.s0(indc(k),:)+post.profil0d.s0m(indc(k),:));
		end
		%ff3c = min(ticxs.fluxC(k,:)) ./ min(fluxc);
		%fluxc = ff3c .* fluxc;
		semilogy(Rti, fluxc,'b',...
		      ticxs.pos(k,:),ticxs.fluxC(k,:),'or');
		if k == 1
		    legend('Metis','Ticxs');
		end
		hold on
		Rti = post.profil0d.Raxe(indc(k),:) - amat(indc(k)) .* post.profil0d.xli;
		semilogy(Rti, fluxc ,'b');
		% creation de vecteur pour le plot des erreur
        	tion=[ticxs.fluxC(k,:);ticxs.fluxC(k,:);ones(size(ticxs.fluxC(k,:)))*NaN];
        	tion=tion(:);
        	pos=[ticxs.pos(k,:)+ticxs.erpos(k,:);ticxs.pos(k,:)-ticxs.erpos(k,:);ones(size(ticxs.pos(k,:)))*NaN];
        	pos=pos(:);
		semilogy(pos,tion,'r');
        	tion=[ticxs.fluxC(k,:)+ticxs.erfluxC(k,:);ticxs.fluxC(k,:)-ticxs.erfluxC(k,:);ones(size(ticxs.fluxC(k,:)))*NaN];
        	tion=tion(:);
        	pos=[ticxs.pos(k,:);ticxs.pos(k,:);ones(size(ticxs.pos(k,:)))*NaN];
        	pos=pos(:);
		semilogy(pos,tion,'r');
		
		if ligne == m
		    xlabel('R (m)');
		end
		if col == 1
		    ylabel('Flux C (eV)');
		end
		if k == 1
		    title(sprintf('Tore Supra, shot # %d @ %g s',fix(shot),post.profil0d.temps(indc(k))));
		else
		    title(sprintf('@ %g s',post.profil0d.temps(indc(k))));
		end
		set(gca,'ButtonDownFcn','zdataplot(''extrait'');');
		
		col = col + 1;
		if col > n
		  col = 1;
		  ligne = ligne +1;
		end
		      
	
		      
        end
        drawnow


	

end
    
