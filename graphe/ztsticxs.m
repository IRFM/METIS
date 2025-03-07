% script de plot de la temperature ionique
% New version 10/01/08 : bug corrected on the mapping of data.prof.ti on the equi.R grid (FI)

shot  = param.from.shot.num;

%ticxs = cgcgettrait(shot,'ticxs');
%if isempty(ticxs.times)
[Ti,tTi]=tsbase(shot,'GTICXS');
[rTi,tTi]=tsbase(shot,'SRAYCXS');
[Tiprof,tprof]=tsbase(shot,'gproftifit_');
rprof = linspace(0,1,21);
if isempty(Ti)
  ticxs = load(sprintf('/home/sccp/gcbd/fenzi/CXRS/ks4fit_java/analyses/TICXS_%d',fix(shot)));
else
  ticxs.pos = rTi(:,1);
  [erTi,tTi]=tsbase(shot,'SERRAYCXS');
  [eTi,tTi]=tsbase(shot,'GERTICXS');
  ticxs.erpos = erTi(:,1);
  ticxs.tcxs = tTi;
  ticxs.Ticxs = Ti;
  ticxs.erTicxs = eTi;
end
%end

if ~isempty(ticxs)	
	
	
	% recherche des indices dans la simulation cronos
	indl = (1:length(data.gene.temps))';
	indc = interp1(data.gene.temps,indl,ticxs.tcxs,'nearest');
	% creation de la figure (une par temps)	
	for k = 1:length(indc)
		h = findobj(0,'type','figure','tag',sprintf('zticxs%d',k));
		if isempty(h)
			h=figure('tag',sprintf('zticxs%d',k));
		else
			figure(h);
		end
		clf
		set(h,'defaultaxesfontsize',16,'defaultaxesfontweight', ...
		      'bold','defaultaxesfontname','times', ...
                      'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
		
		% changement de coordonnees      
		Rti = squeeze(data.equi.R(indc(k),:,1));
		if all(~isfinite(Rti))
			Rti = (max(data.geo.R(indc(k),:)) + min(data.geo.R(indc(k),:))) ./2  + ...
			      (max(data.geo.R(indc(k),:)) -  min(data.geo.R(indc(k),:))) ./2 .* ...
			      param.gene.x + 0.1 .* (1 - param.gene.x .^ 2);
		        disp('no equilibirum for this time');
		end
		ti_equigrid = interp1(param.gene.x.*data.equi.rhomax(indc(k)),data.prof.ti(indc(k),:),data.equi.rhoRZ(indc(k),:));
		plot(Rti,ti_equigrid,'b',ticxs.pos,ticxs.Ticxs,'or');
		legend('Cronos','Ticxs');
		hold on
		indhfs  = floor(size(data.equi.R,3)./2);
		Rti = squeeze(data.equi.R(indc(k),:,indhfs));
		if all(~isfinite(Rti))
			Rti = (max(data.geo.R(indc(k),:)) + min(data.geo.R(indc(k),:))) ./2  - ...
			      (max(data.geo.R(indc(k),:)) -  min(data.geo.R(indc(k),:))) ./2 .* ...
			      param.gene.x + 0.1 .* (1 - param.gene.x .^ 2);
		        disp('no equilibrium for this time');
		end
		plot(Rti,ti_equigrid,'b');
		
		% creation de vecteur pour le plot des erreur
        	tion=[ticxs.Ticxs';ticxs.Ticxs';ones(size(ticxs.Ticxs'))*NaN];
        	tion=tion(:);
            if size(ticxs.pos,1) == 1
        	  pos=[ticxs.pos'+ticxs.erpos';ticxs.pos'-ticxs.erpos';ones(size(ticxs.pos'))*NaN];
        	  pos=pos(:);
            else
        	  pos=[ticxs.pos+ticxs.erpos;ticxs.pos-ticxs.erpos;ones(size(ticxs.pos))*NaN];
        	  pos=pos(:);                
            end
            long = length(ticxs.Ticxs);
        for i=1:long   
		  plot(pos([i i+long]),tion([i i+long]),'r');
        end
        	tion=[ticxs.Ticxs'+ticxs.erTicxs';ticxs.Ticxs'-ticxs.erTicxs';ones(size(ticxs.Ticxs'))*NaN];
        	tion=tion(:);
            if size(ticxs.pos,1) == 1
            	pos=[ticxs.pos';ticxs.pos';ones(size(ticxs.pos'))*NaN];
        	    pos=pos(:);
            else
             	pos=[ticxs.pos;ticxs.pos;ones(size(ticxs.pos))*NaN];
        	    pos=pos(:);
               
            end
    for i=1:long   
			plot(pos([i i+long]),tion([i i+long]),'b');
       end
 		
		xlabel('R (m)');
		ylabel('T_i (eV)');
		title(sprintf('Tore Supra, shot # %d @ %g s',fix(shot),data.gene.temps(indc(k))));
		      
        end

%
% avec tprof
%
	
if ~isempty(Tiprof)

	h = findobj(0,'type','figure','tag',sprintf('zticxs_tprof%d',k));
	if isempty(h)
		h=figure('tag',sprintf('zticxs_tprof%d',k));
	else
		figure(h);
	end
	clf
	set(h,'defaultaxesfontsize',16,'defaultaxesfontweight', ...
		      'bold','defaultaxesfontname','times', ...
                      'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])

	for k = 1:length(indc)
	  indprof = iround(tprof,data.gene.temps(indc(k)));
 	  plot(param.gene.x,data.prof.ti(indc(k),:),'b',rprof,Tiprof(indprof,:)*1e3,'or');
	  legend('Cronos','Tprof');
        end



end

end
    