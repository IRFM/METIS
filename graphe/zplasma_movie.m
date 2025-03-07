% ZPLASMA_MOVIE  cree un film a partir des donnees cronos
%------------------------------------------------------------
% fichier zplama_movie.m ->  zplasma_movie
%
% fonction Matlab 5 :
% Cette fonction cree un film a partir des donnees cronos.
% Il s'agit d'un exmple qui doit etre modifie en fonction de ce que l'on veut tracer
%
% syntaxe  :
%
% pour cree un fichier avi :
%
%	filename = zplasma_movie(param,data);
%       + donnee le nom du fichier dans l'interface
%
% pour cree un 'movie' au sens matlab :
%	mov = zplasma_movie(param,data);
%       + ne pas mettre de nom de fichier dans l'interface
%
% pour visualise le film sans le cree :
%	zplasma_movie(param,data);
%
%
% entrees :
%	param = structure de parametres cronos
%	data  = structure de donnees cronos
%
% sorties : 
%	filename = nom du fichier cree
%       mov      = objet matlab movie
%
%
% fonction ecrite par J-F Artaud
% version vcvs,  (cf. detail sur le serveur cvs).
%
%
% liste des modifications : cf. serveur cvc repository : gcgc, module : Projet_Cronos
%
%--------------------------------------------------------------
%
function filename=zplasma_movie(param,data)

% mov
mov =[];

% nom du fichier AVI
if nargout > 0
	[filename, pathname] = uiputfile('*.avi', 'Name of AVI file');
else
	filename = 0;
	pathname = 0;
end
if ischar(filename)
	[voidpath,filename] = fileparts(filename);
	filename = fullfile(pathname,strcat(filename,'.avi'))
	mov = avifile(filename,'compression','none');
end

% figure
hf = findobj(0,'type','figure','tag','zplasma_movie');
if isempty(hf)
	hf = figure;
else
	figure(hf);
end
clf
set(hf,'tag','zplasma_movie','DoubleBuffer','on','color',[0 0 0], ...
            'defaultaxesfontsize',16,'defaultaxesfontweight','bold', ...
	    'defaultaxesfontname','times','defaulttextcolor',[1 1 1],'defaultaxescolor',[0 0 0], ...
	    'defaultaxeszcolor',[1 1 1],'defaultaxesxcolor',[1 1 1],'defaultaxesycolor',[1 1 1]);

title('Set the rigth figure size and, after,  strike a key');
[xv,yv] = ginput(1);
title('Now, do not resize the figure !');
pause(1)
    
hc=uicontrol(hf,'style','radio','tag','onpause','value',0,'string','pause');


% pas
pas = max(1,round(length(param.gene.x) / 21));

% maximum
if param.gene.k > param.gene.kmin
	indok =param.gene.kmin:param.gene.k;
else
	indok =param.gene.kmin:param.gene.kmax;
end
qmax = max(data.prof.q(indok,end));
pmax = min(1e9,max(data.prof.ptot(indok,1)));
tmax = min(100e3,max(max(data.prof.te(indok,1) ,data.prof.ti(indok,1))));
nmax = min(1e22,max(data.prof.ne(indok,1)));
jmax = min(1e9,max(max(abs(data.prof.jmoy(indok,:)))));
negr   = 1e20 .* (data.gene.ip /1e6) ./ (pi.* data.equi.a(:,end) .^ 2);
pfus = data.source.fus.el(indok,:) + data.source.fus.ion(indok,:);
fmax  = max(max(pfus));

% (R,Z)
R = double(data.equi.R(indok,:,:));
Z = double(data.equi.Z(indok,:,:));

% echelle
Rmin = min(R(:)) .* 0.95;
Rmax = max(R(:)) .* 1.05;
Zmin = min(Z(:)) .* 1.05;
Zmax = max(Z(:)) .* 1.05;

% nombre de frame
nbf = length(indok);


% boucle sur les temps
first = 1;
for k =1:nbf
  % temps courant
  tc = data.gene.temps(k);
  colmem =[];
  figure(hf)
  hold off
  cla
  set(gca,'box','on');
  colormap('gray')
  hold on
  % le fond  (mise en evidence du rayonnement)
  pfond = data.source.totale.el(k,length(param.gene.x):-pas:2)' + ...
	  data.source.totale.ion(k,length(param.gene.x):-pas:2)';
  pfond = - pfond .* (pfond < 0);
  Rfond = squeeze(double(data.equi.R(k,length(param.gene.x):-pas:2,:)));
  Zfond = squeeze(double(data.equi.Z(k,length(param.gene.x):-pas:2,:)));
  pfond = pfond * ones(1,size(Rfond,2));
  hsurf = surf(Rfond',Zfond',0.*Zfond',pfond','linestyle','none', ...
               'facecolor','interp','edgecolor','none');
  


  
  % boucle sur les surface de flux
  for l = length(param.gene.x):-pas:2
  	if any((param.compo.a == 3)&(param.compo.z == 1))
		rouge      = max(0,min(1,(max(data.prof.te(k,l),data.prof.ti(k,l)) ./ tmax)));
		bleu       = 1 - max(0,min(1,(pfus(k,l)) ./ fmax));
		vert       = 1 - max(0,tanh((data.prof.ne(k,l) - data.prof.ne(k,end))./ nmax));
		couleur = [rouge vert bleu] ;
		% epais = max(0.5,6.5 - profli.qjli(k,l));
		if data.prof.q(k,l) <= 1
			epais = 6;
		elseif data.prof.q(k,l) <= 1.5
			epais = 5;
		elseif data.prof.q(k,l) <= 2
			epais = 4;
		elseif data.prof.q(k,l) <= 2.5
			epais = 3;
		elseif data.prof.q(k,l) <= 3
			epais = 2;
		elseif data.prof.q(k,l) <= 5
			epais = 1;
		else
			epais = 0.5;
		end

	else
		rouge      = max(0,min(1,(max(data.prof.te(k,l),data.prof.ti(k,l)) ./ tmax)));
		bleu       = 1 - max(0,min(1,data.prof.ptot(k,l) ./ pmax));
		vert       = 1 - max(0,tanh((data.prof.ne(k,l) - data.prof.ne(k,end))./ nmax));
		couleur = [rouge vert bleu] ;
		epais = max(0.5,6.5 - data.prof.q(k,l));
	end
		
	plot(squeeze(double(data.equi.R(k,l,:))),squeeze(double(data.equi.Z(k,l,:))), ...
	    'color',couleur,'linewidth',epais);
	if  data.prof.q(k,l) <= 1
		plot(squeeze(double(data.equi.R(k,l,1:3:end))),squeeze(double(data.equi.Z(k,l,1:3:end))), ...
		'linestyle',':','color',1 - couleur,'linewidth',epais);
	end
	if isempty(colmem)
		colmem = couleur;
	end
  end
  
  if data.geo.mode == 2
  	Rs = double(data.geo.R(k,:));
	Zs = double(data.geo.Z(k,:));
  	plot(Rs,Zs,'color',1-colmem,'linewidth',0.5,'linestyle','-.');	
  else
  	plot(squeeze(double(data.equi.R(k,end,:))),squeeze(double(data.equi.Z(k,end,:))), ...
	     'color',1-colmem,'linewidth',0.5,'linestyle','-.')
  end
  plot(data.equi.raxe(k,end),mean(squeeze(double(data.equi.Z(k,2,:)))),'marker','+','linestyle','none','color',1-couleur);
  
  xlabel('R (m)')
  ylabel('Z (m)')
  if k == 1
   	set(gca,'xlim',[Rmin,Rmax],'ylim',[Zmin,Zmax]);
 	axis('equal');
	drawnow
	xl = get(gca,'xlim');
	yl = get(gca,'ylim');
	tl = get(gca,'TickLength');
  else
  	set(gca,'xlim',xl,'ylim',yl); 
  end


  title(sprintf('time = %g s [%s @ %g]',data.gene.temps(k),param.from.machine,param.from.shot.num));
  ht = text(0,0,'123456789012345678');
  ex = get(ht,'extent');
  delete(ht);
  
  Rmin = xl(1)+tl(1).* abs(diff(xlim));
  Rmax = xl(2)-tl(1).* abs(diff(xlim)); 
  Zmin = yl(1)+tl(2).* abs(diff(ylim));
  Zmax = yl(2)-tl(2).* abs(diff(ylim));

  text(Rmin,Zmax - 1.* ex(4),sprintf(' I_p = %6.3g MA',data.gene.ip(k)./1e6));
  text(Rmin,Zmax - 2 .* ex(4),sprintf(' I_b_o_o_t = %6.3g MA',data.gene.iboot(k)./1e6));
  text(Rmin,Zmax - 3 .* ex(4),sprintf(' I_c_d = %6.3g MA',data.gene.icd(k)./1e6));
  text(Rmin,Zmax - 4 .* ex(4),sprintf(' I_n_i = %6.3g MA',data.gene.ini(k)./1e6));
  text(Rmin,Zmax - 5 .* ex(4),sprintf(' V_l_o_o_p = %6.3g V',data.gene.vloop(k)));
  
  text(Rmin,Zmin + 7.* ex(4),sprintf(' P_N_B_I = %6.3g MW',data.gene.paddidn(k)./1e6));
  text(Rmin,Zmin + 6 .* ex(4),sprintf(' P_I_C_R_H = %6.3g MW',data.gene.paddfci(k)./1e6));
  text(Rmin,Zmin + 5 .* ex(4),sprintf(' P_L_H = %6.3g MW',data.gene.paddhyb(k)./1e6));
  text(Rmin,Zmin + 4 .* ex(4),sprintf(' P_E_C_R_H = %6.3g MW',data.gene.paddfce(k)./1e6));
  text(Rmin,Zmin + 3 .* ex(4),sprintf(' P_f_u_s = %6.3g MW',data.gene.paddfus(k)./1e6));
  text(Rmin,Zmin + 2 .* ex(4),sprintf(' P_l_i_n_e = %6.3g MW',data.gene.prad(k)./1e6));
  text(Rmin,Zmin + 1 .* ex(4),sprintf(' P_b_r_e_m = %6.3g MW',data.gene.pbrem(k)./1e6));
  text(Rmin,Zmin + 0 .* ex(4),sprintf(' P_s_y_n_c = %6.3g MW',data.gene.pcyclo(k)./1e6));

  text(Rmax - ex(3), Zmax - 1 .* ex(4),sprintf('T_e_0 = %6.3g keV',data.prof.te(k,1)./1e3));
  text(Rmax - ex(3), Zmax - 2 .* ex(4),sprintf('T_i_0 = %6.3g keV',data.prof.ti(k,1)./1e3));
  text(Rmax - ex(3), Zmax - 3 .* ex(4),sprintf('n_e_0 = %6.3g 1e19 m^-^3',data.prof.ne(k,1)./1e19));
  text(Rmax - ex(3), Zmax - 4 .* ex(4),sprintf('n_b_a_r = %6.3g 1e19 m^-^3',data.gene.nbar(k)./1e19));
  text(Rmax - ex(3), Zmax - 5 .* ex(4),sprintf('f_g_r = %6.3g ',data.gene.nbar(k)./negr(k)));

  fex = 0.7;
  text(Rmax - ex(3).*fex, Zmin + 5 .* ex(4),sprintf('beta_P = %6.3g ',data.gene.betap(k)));
  text(Rmax - ex(3).*fex, Zmin + 4 .* ex(4),sprintf('l_i = %6.3g ',data.gene.li(k)));
  text(Rmax - ex(3).*fex, Zmin + 3 .* ex(4),sprintf('q_a = %6.3g ',data.prof.q(k,end)));
  text(Rmax - ex(3).*fex, Zmin + 2 .* ex(4),sprintf('q_0 = %6.3g ',data.prof.q(k,2)));
  text(Rmax - ex(3).*fex, Zmin + 1 .* ex(4),sprintf('q_m_i_n = %6.3g ',min(data.prof.q(k,:))));
  indqmin = max(find(data.prof.q(k,:)== min(data.prof.q(k,:))));
  text(Rmax - ex(3).*fex, Zmin + 0 .* ex(4),sprintf('x_m_i_n = %6.3g ',param.gene.x(indqmin)));

  drawnow
  
  % ajout d'une image
  if ischar(filename)
  	 F = getframe(hf);
         mov = addframe(mov,F);
  elseif nargout > 0
   	 F = getframe(hf);
 	 mov(k) =F;
  elseif (get(hc,'value') == 1)
  	set(gca,'zcolor',[0 0 0],'xcolor',[0 0 0],'ycolor',[0 0 0],'color',[1 1 1]);
	set(get(gca,'xlabel'),'color',[0 0 0]);
	set(get(gca,'ylabel'),'color',[0 0 0]);
	set(get(gca,'title'),'color',[0 0 0]);
	set(findobj(gcf,'type','text'),'color',[0 0 0]);
	set(gcf,'color',[1 1 1]);
	map = colormap;
	colormap(1-map);
	
  	waitfor(hc,'value');
	
   	set(gca,'zcolor',[1 1 1],'xcolor',[1 1 1],'ycolor',[1 1 1],'color',[0 0 0]);
	set(get(gca,'xlabel'),'color',[1 1 1]);
	set(get(gca,'ylabel'),'color',[1 1 1]);
	set(get(gca,'title'),'color',[1 1 1]);
	set(findobj(gcf,'type','text'),'color',[1 1 1]);
	set(gcf,'color',[0 0 0]);
	colormap('gray');
  else
	 pause(30/nbf);
  end

end

% fermeture du fichier
if ischar(filename)
	mov = close(mov);
else
	filename = mov;
end

pause(1)
set(gca,'zcolor',[0 0 0],'xcolor',[0 0 0],'ycolor',[0 0 0],'color',[1 1 1]);
set(findobj(gcf,'type','text'),'color',[0 0 0]);
set(gcf,'color',[1 1 1]);
map = colormap;
colormap(1-map);
