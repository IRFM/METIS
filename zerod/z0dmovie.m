% Z0DMOVIE  cree un film a partir des donnees de Metis
%------------------------------------------------------------
% fichier z0dmovie.m ->  z0dmovie
%
% fonction Matlab 5 :
% Cette fonction cree un film a partir des donnees cronos.
% Il s'agit d'un exmple qui doit etre modifie en fonction de ce que l'on veut tracer
%
% syntaxe  :
%
% pour cree un fichier avi :
%
%	filename = z0dmovie(post);
%       + donnee le nom du fichier dans l'interface
%
% pour cree un 'movie' au sens matlab :
%	mov = z0dmovie(post);
%       + ne pas mettre de nom de fichier dans l'interface
%
% pour visualise le film sans le cree :
%	zplasma_movie(post);
%
%
% entrees :
%	post  = structure de donnees post de CRONOS/Metis
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
% remarque : compression sous linux =  ! mencoder input.avi -o output_comp.avi -ovc lavc
%
%--------------------------------------------------------------
%
function filename=z0dmovie(post,temps)

persistent open_gl_software
if isempty(open_gl_software) 
  if exist('opengl','builtin')
      void = opengl('DATA');
      open_gl_software = void.Software;
  else
    open_gl_software = 0;
  end
end

if isappdata(0,'GUI_GLOBAL_BASIC_MODE')
   basic_mode = getappdata(0,'GUI_GLOBAL_BASIC_MODE');
else
   basic_mode = 0;
end
% mov
mov =[];

% nom du fichier AVI
if nargout > 0
	[filename, pathname] = uiputfile('*.avi', 'Name of AVI file');
    if isequal(filename,0) || isequal(pathname,0)
            return
    end
else
	filename = 0;
	pathname = 0;
end
if ischar(filename)
	[voidpath,filename] = fileparts(filename);
	filename = fullfile(pathname,strcat(filename,'.avi'))
	if exist('VideoWriter')
	  vidObj = VideoWriter(filename);
	  open(vidObj);
	else
	  mov = avifile(filename,'compression','none');
	end
end

if ~isfield(post.z0dinput.option,'te_max')
    te_max = 100;
else
    te_max = post.z0dinput.option.te_max;
end


% create legend figure
% figure
hf = findobj(0,'type','figure','tag','z0dmovie_legend');
if isempty(hf)
	hf = figure;
else
	figure(hf);
end
clf
pos = get(hf,'position');
set(hf,'Color',[1 1 1],'defaultaxesFontSize',18,'position',[1 1 pos(3) pos(4)],'defaultaxescolor',[0 0 0],'tag','z0dmovie_legend');

for k=0:0.1:3
    betap = k;
    factcol = tanh(betap);
    rouge   = factcol; 
    bleu    = 1 - rouge;
    vert    = sqrt(1 - bleu .^ 2 - rouge .^ 2);
    couleur = [rouge vert bleu] ;
    for l=1:7 
      switch l
      case 1
 	epais = 6;
 	q     = [0.5,1];
      case 2
  	epais = 5;
 	q     = [1,1.5];
      case 3
   	epais = 4;
 	q     = [1.5,2];
      case 4
   	epais = 3;
 	q     = [2,2.5];
      case 5
   	epais = 2;
 	q     = [2.5,3];
      case 6
    	epais = 1;
 	q     = [3,5];
      case 7
    	epais = 0.5;
 	q     = [5,8];
      end
      if l > 1
	q = [q(1) + 0.1,q(2)];
      end
      if l < 7
	q = [q(1),q(2) - 0.1];      
      end
      plot(q,[betap,betap],'linewidth',epais,'color',couleur); 
      hold on
      if  all(q <= 1)
	  plot(q,[betap,betap],'linestyle',':','color',1 - couleur);
      end
    end
end
xlabel('safety factor interval');
ylabel('local \beta_P');  
if (basic_mode == 0) && (verLessThan('matlab','8.6') || (open_gl_software == 0))
  title(sprintf('video''s legend\n(black-gray-white background represent line radiation intensity)'));
else
  title('video''s legend');
end
axis([0 8.5 -0.1 3.1]);



% figure
hf = findobj(0,'type','figure','tag','z0dmovie');
if isempty(hf)
	hf = figure;
else
	figure(hf);
end
clf
set(hf,'tag','z0dmovie','DoubleBuffer','on','color',[0 0 0], ...
            'defaultaxesfontsize',16,'defaultaxesfontweight','bold', ...
	    'defaultaxesfontname','times','defaulttextcolor',[1 1 1],'defaultaxescolor',[0 0 0], ...
	    'defaultaxeszcolor',[1 1 1],'defaultaxesxcolor',[1 1 1],'defaultaxesycolor',[1 1 1],'resize','on');


title('Set the rigth figure size and, after,  strike a key');
[xv,yv] = ginput(1);
set(hf,'resize','off');
title('Now, movie will start');
pause(1)


% pas
profli = post.profil0d;
zs     = post.zerod;
geo    = post.z0dinput.geo;
if ~isfield(geo,'z0')
	disp('undefined z0 : z0 set to 0');
	geo.z0 = 0.*geo.a;
end
pas = max(1,round(length(profli.xli) / 10));

% (R,Z) separatrice
if isfield(post.z0dinput.exp0d,'Rsepa') & isfield(post.z0dinput.exp0d,'Zsepa')
	Rs = post.z0dinput.exp0d.Rsepa;
	Zs = post.z0dinput.exp0d.Zsepa +  geo.z0 * ones(1,size(Rs,2));
        maskrmax  = (Rs == (max(Rs,[],2) * ones(1,size(Rs,2))));
        z0   = sum(Zs .* maskrmax,2) ./ sum(maskrmax,2);
else

	td  = asin(max(0,min(1,geo.d)));
	u   = linspace(0,2.*pi,201);
	vu  = ones(size(u));
	vt  = ones(size(td));
	Rs  = geo.R *vu + (geo.a * vu) .* cos(vt * u + td * sin(u));
	Zs  = (geo.a .* geo.K) * sin(u) + geo.z0 * ones(1,size(Rs,2));
	z0  = geo.z0;

end


% limiter
rwall = [];
zwall = [];
if isfield(post.z0dinput.option,'first_wall')
    [pfw,ffw] = fileparts(post.z0dinput.option.first_wall);
    filename_wall = fullfile(pfw,sprintf('%s.mat',ffw));
    if exist(filename_wall,'file')
            wall = load(filename_wall);
            if isfield(wall,'R') && isfield(wall,'Z')
                rwall = wall.R(:);
                zwall = wall.Z(:);
            end
    end
end
if isempty(rwall)
    switch post.z0dinput.machine
    case 'ITER'
        if  exist('iterwall')
            [rwall,zwall] = iterwall;
        else
            rwall = [];
            zwall = [];
        end    
    otherwise
        if ~isempty(strfind(upper(post.z0dinput.machine),'WEST')) 
		[rwall,zwall] = west_limiter(post.z0dinput.shot);
		try
		    walldata = Get_Paroi_WEST(post.z0dinput.shot,post.z0dinput.cons.temps);
		    rwall    = walldata.Rparoi;
		    zwall    = walldata.Zparoi;
		    rwall(:,end+1) = rwall(:,1);
		    zwall(:,end+1) = zwall(:,1);                  
		end

        elseif ~isempty(strfind(post.z0dinput.machine,'JT-60SA')) || ~isempty(strfind(post.z0dinput.machine,'JT60-SA'))
                wall = load('jt60sawall');
                rwall = wall.wall.data(1:end,1);
                zwall = wall.wall.data(1:end,2);
        else
            rwall = [];
            zwall = [];
        end
    end
end

% maximum
qmax = max(profli.qjli(:,end));
pmax = min(1e9,max(profli.ptot(:,1)));
tmax = min(te_max,max(max(profli.tep(:,1) ,profli.tip(:,1))));
nmax = min(1e22,max(profli.nep(:,1))) - max(0,min(profli.nep(:,1)));
jmax = min(1e9,max(max(abs(profli.jli(:,:)))));
negr   = 1e20 .* (zs.ip /1e6) ./ (pi.* geo.a .^ 2);
fmax  = max(max(profli.pfus));


% echelle
Rmin = min(Rs(:)) .* 0.95;
Rmax = max(Rs(:)) .* 1.05;
Zmin = min(Zs(:)) .* 1.05;
Zmax = max(Zs(:)) .* 1.05;

% prise en compte du mur
if ~isempty(rwall) && ~isempty(zwall)
  Rmin = min(Rmin ,min(rwall(:)) .* 0.95);
  Rmax = max(Rmax, max(rwall(:)) .* 1.05);
  Zmin = min(Zmin, min(zwall(:)) .* 1.05);
  Zmax = max(Zmax, max(zwall(:)) .* 1.05);
end

% si 2 entrees
if nargin < 2
	% nombre de frames
	nbf = length(profli.temps);
	nbini = 1;
	temps = profli.temps;
else
	nbini = find(profli.temps >= min(temps),1);
	nbf   = max(find(profli.temps <= max(temps)));
        nbini = min(nbf,nbini);
end
llist = [2,4,6,8,10,12,14,16,18,20,21];
u    = linspace(0,2.*pi,65);

% bouton de controle
if length(temps)>1
	hc=uicontrol(hf,'style','radio','tag','onpause','value',0,'string','pause');
else
	hc=uicontrol(hf,'style','radio','tag','onpause','value',0,'string','pause','visible','off');
end

% boucle sur les temps
first = 1;
for kf =nbini:nbf
    % temps courant
    tc = profli.temps(kf);
    k  = min(find(tc <= post.zerod.temps));
    if isempty(k)
        k= length(post.zerod.temps);
    end
    
    figure(hf)
    hold off
    cla
    colormap('gray')
    set(gca,'box','on');
    hold on
    
    % boucle sur les surface de flux
    %clear Rfond Zfond
    for ll = 1:length(llist)
        l =llist(ll);
        switch post.z0dinput.option.gaz
            case 3
                %  		rouge      = max(0,min(1,(max(profli.tep(kf,l),profli.tip(kf,l)) ./ tmax)));
                %  		bleu       = 1 - max(0,min(1,(profli.pfus(kf,l)) ./ fmax));
                %                  vert       = 1 - max(0,tanh(profli.nep(kf,l) ./ max(negr)));
                %  		couleur = [rouge vert bleu] ;
                factcol = tanh(2 .* (4.*pi .* 1e-7) .* profli.ptot(kf,l) ./ profli.bpol(kf,end) .^ 2);
                rouge   = factcol;
                bleu    = 1 - rouge;
                vert    = sqrt(1 - bleu .^ 2 - rouge .^ 2);
                couleur = [rouge vert bleu] ;
                % epais = max(0.5,6.5 - profli.qjli(kf,l));
                if profli.qjli(kf,l) <= 1
                    epais = 6;
                elseif profli.qjli(kf,l) <= 1.5
                    epais = 5;
                elseif profli.qjli(kf,l) <= 2
                    epais = 4;
                elseif profli.qjli(kf,l) <= 2.5
                    epais = 3;
                elseif profli.qjli(kf,l) <= 3
                    epais = 2;
                elseif profli.qjli(kf,l) <= 5
                    epais = 1;
                else
                    epais = 0.5;
                end
                
            otherwise
                %  		rouge      = max(0,min(1,(max(profli.tep(kf,l),profli.tip(kf,l)) ./ tmax)));
                %  		bleu       = 1 -  max(0,min(1,profli.ptot(kf,l) ./ pmax));
                %  		% vert       = 1 - max(0,tanh((profli.nep(kf,l) - profli.nep(kf,end))./ nmax));
                %  		vert       = 1 - max(0,tanh(profli.nep(kf,l) ./ max(negr)));
                %  		couleur = [rouge vert bleu] ;
                factcol = tanh(2 .* (4.*pi .* 1e-7) .* profli.ptot(kf,l) ./ profli.bpol(kf,end) .^ 2);
                rouge   = factcol;
                bleu    = 1 - rouge;
                vert    = sqrt(1 - bleu .^ 2 - rouge .^ 2);
                couleur = [rouge vert bleu] ;
                
                if profli.qjli(kf,l) <= 1
                    epais = 6;
                elseif profli.qjli(kf,l) <= 1.5
                    epais = 5;
                elseif profli.qjli(kf,l) <= 2
                    epais = 4;
                elseif profli.qjli(kf,l) <= 2.5
                    epais = 3;
                elseif profli.qjli(kf,l) <= 3
                    epais = 2;
                elseif profli.qjli(kf,l) <= 5
                    epais = 1;
                else
                    epais = 0.5;
                end
        end
        %epais = max(0.5,6.5 - profli.qjli(kf,l));
        
        td  = asin(max(0,min(1,geo.d(k) .* (profli.xli(l) .^ 3))));
        R   = profli.Raxe(kf,l) + geo.a(k) .* profli.xli(l).* cos(u + td * sin(u));
        Z   = geo.a(k) .* profli.xli(l) .* profli.kx(kf,l) * sin(u)+ z0(k) ;
        
        % deformation
        cx   = (R - profli.Raxe(kf,l)) + sqrt(-1) .* (Z - z0(k));
        thx  = unwrap(angle(cx),[],2);
        rhox = abs(cx);
        thx(thx<0) = thx(thx<0) + 2 .* pi;
        [thx,indx] = sort(thx);
        rhox       = rhox(indx);
        
        % separtrice a ce temps version analytique
        tdl  = asin(max(0,min(1,geo.d(k))));
        Rl   = profli.Raxe(kf,end) + geo.a(k) .* cos(u + tdl * sin(u));
        Zl   = geo.a(k).* profli.kx(kf,end) * sin(u)+ z0(k);
        cl   = (Rl - profli.Raxe(kf,l)) + sqrt(-1) .* (Zl - z0(k));
        thl  = unwrap(angle(cl),[],2);
        rhol = abs(cl);
        thl(thl<0) = thx(thl<0) + 2 .* pi;
        [thl,indl] = sort(thl);
        rhol       = rhol(indl);
        rhol = cat(2,rhol,rhol,rhol);
        thl = cat(2,thl -2.*pi,thl,thl+2.*pi);
        indnok = find(diff(thl)<=0);
        thl(indnok) =[];
        rhol(indnok)  = [];
        
        % separtrice a ce temps complete
        cc   = (Rs(k,:) - profli.Raxe(kf,l)) + sqrt(-1) .* (Zs(k,:) - z0(k));
        thc  = unwrap(angle(cc),[],2);
        rhoc = abs(cc);
        thc(thc<0) = thc(thc<0) + 2 .* pi;
        [thc,indc] = sort(thc);
        rhoc       = rhoc(indc);
        rhoc = cat(2,rhoc,rhoc,rhoc);
        thc = cat(2,thc -2.*pi,thc,thc+2.*pi);
        indnok = find(diff(thc)<=0);
        thc(indnok) =[];
        rhoc(indnok)  = [];
        
        
        % regle de trois, mise a l'echelle
        if post.z0dinput.option.morphing > 0
            morf    = (1 - profli.xli(l).^ post.z0dinput.option.morphing) +  ...
                profli.xli(l) .^ post.z0dinput.option.morphing .* pchip(thc',rhoc',thx')'  ./ pchip(thl',rhol',thx')';
        else
            morf = ones(size(thx));
        end
        rhox = rhox .* morf;
        Rx   = rhox .* cos(thx) + profli.Raxe(kf,l);
        Zx   = rhox .* sin(thx)+ z0(k);
        Rx(end+1) = Rx(1);
        Zx(end+1) = Zx(1);
        
        % tracer de la surface de flux
        plot(Rx,Zx,'color',couleur,'linewidth',epais);
        hold on
        
        if  profli.qjli(kf,l) <= 1
            plot(Rx,Zx,'linestyle',':','color',1 - couleur);
        end
        
        Rfond(l,:) = Rx;
        Zfond(l,:) = Zx;
        
        if l == 2
            couleur_x = 1 - couleur;
        end
        
    end
    plot(Rs(k,:),Zs(k,:),'color',[0 0 0],'linewidth',0.5,'linestyle','-.')
    plot(geo.R(k),z0(k),'+','color',couleur_x);
    
    if ~isempty(rwall) && ~isempty(zwall)
        if all(size(rwall) > 1) && (size(rwall,1) == length(temps))
            plot(rwall(k,:),zwall(k,:),'m')
        else
            plot(rwall,zwall,'m')
        end
    end
    
    
    % le fond  (mise en evidence du rayonnement)
    if basic_mode == 0
        %if verLessThan('matlab','8.6') || (open_gl_software == 0)
        Rfond = Rfond(llist,:);
        Zfond = Zfond(llist,:);
        pfond = post.profil0d.prad(kf,llist)' + post.profil0d.pioniz(kf,llist)';
        pfond = pfond .* (pfond > 0);
        pfond = pfond * ones(1,size(Rfond,2));
        hsurf = surf(Rfond',Zfond',0.*Zfond',pfond','linestyle','none', ...
            'facecolor','interp','edgecolor','none');
        if ~verLessThan('matlab','8.6')
            alpha( hsurf,.33);
        end
        %end
    end
    
    xlabel('R (m)')
    ylabel('Z (m)')
    if first == 1
        ok =0;
        rrr = 2;
        while (ok ==0)
            Rm = (Rmin + Rmax) ./ 2;
            dRm = Rmax - Rmin;
            Zm = (Zmin + Zmax) ./ 2;
            dZm = Zmax - Zmin;
            ddm = max(dRm,dZm)./rrr;
            set(gca,'xlim',Rm + [-ddm,ddm],'ylim',Zm + [-ddm,ddm]);
            axis('equal');
            drawnow
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            tl = get(gca,'TickLength');
            if (min(xl) >= Rmin) | (max(xl) <= Rmax) |(min(yl) >= Zmin) | (max(yl)<= Zmax)
                rrr = rrr .* 0.9 ;
            else
                ok = 1;
            end
        end
    else
        set(gca,'xlim',xl,'ylim',yl);
    end
    
    title(sprintf('time = %g s',profli.temps(kf)))
    
    
    ht = text(0,0,'123456789012345678');
    ex = get(ht,'extent');
    delete(ht);
    
    Rmin = xl(1)+tl(1).* abs(diff(xlim));
    Rmax = xl(2)-tl(1).* abs(diff(xlim));
    Zmin = yl(1)+tl(2).* abs(diff(ylim));
    Zmax = yl(2)-tl(2).* abs(diff(ylim));
    
    text(Rmin,Zmax - 1 .* ex(4),sprintf(' I_p = %6.3g MA',zs.ip(k)./1e6));
    text(Rmin,Zmax - 2 .* ex(4),sprintf(' I_b_o_o_t = %6.3g MA',zs.iboot(k)./1e6));
    text(Rmin,Zmax - 3 .* ex(4),sprintf(' I_c_d = %6.3g MA',zs.icd(k)./1e6));
    text(Rmin,Zmax - 4 .* ex(4),sprintf(' I_n_i = %6.3g MA',zs.ini(k)./1e6));
    text(Rmin,Zmax - 5 .* ex(4),sprintf(' V_l_o_o_p = %6.3g V',zs.vloop(k)));
    
    text(Rmin,Zmin + 7.* ex(4),sprintf(' P_N_B_I = %6.3g MW',real(zs.pnbi_th(k))./1e6  .* (real(zs.pnbi_th(k)) > 1) + ...
        imag(zs.pnbi_th(k))./1e6 .* (imag(zs.pnbi_th(k)) > 1)));
    text(Rmin,Zmin + 6 .* ex(4),sprintf(' P_I_C_R_H = %6.3g MW',zs.picrh_th(k)./1e6 .* (zs.picrh_th(k) > 1)));
    if post.z0dinput.option.lhmode == 5
        text(Rmin,Zmin + 5 .* ex(4),sprintf(' P_E_C_R_H_@_1 = %6.3g MW',zs.plh(k)./1e6 .* (zs.plh(k) > 1)));
        text(Rmin,Zmin + 4 .* ex(4),sprintf(' P_E_C_R_H_@_2 = %6.3g MW',zs.pecrh(k)./1e6 .* (zs.pecrh(k) > 1)));
    else
        text(Rmin,Zmin + 5 .* ex(4),sprintf(' P_L_H = %6.3g MW',zs.plh(k)./1e6 .* (zs.plh(k) > 1)));
        text(Rmin,Zmin + 4 .* ex(4),sprintf(' P_E_C_R_H = %6.3g MW',zs.pecrh(k)./1e6 .* (zs.pecrh(k) > 1)));
    end
    text(Rmin,Zmin + 3 .* ex(4),sprintf(' P_{alpha} = %6.3g MW',zs.pfus_th(k)./1e6 .* (zs.pfus_th(k) > 1)));
    text(Rmin,Zmin + 2 .* ex(4),sprintf(' P_l_i_n_e = %6.3g MW',zs.prad(k)./1e6 .* (zs.prad(k) > 1)));
    text(Rmin,Zmin + 1 .* ex(4),sprintf(' P_b_r_e_m = %6.3g MW',zs.pbrem(k)./1e6 .* (zs.pbrem(k) > 1)));
    text(Rmin,Zmin + 0 .* ex(4),sprintf(' P_s_y_n_c = %6.3g MW',zs.pcyclo(k)./1e6 .* (zs.pcyclo(k) > 1)));
    
    text(Rmax - ex(3), Zmax - 1 .* ex(4),sprintf('T_e_0 = %6.3g keV',profli.tep(kf,1)./1e3));
    text(Rmax - ex(3), Zmax - 2 .* ex(4),sprintf('T_i_0 = %6.3g keV',profli.tip(kf,1)./1e3));
    text(Rmax - ex(3), Zmax - 3 .* ex(4),sprintf('n_e_0 = %6.3g 1e19 m^-^3',profli.nep(kf,1)./1e19));
    text(Rmax - ex(3), Zmax - 4 .* ex(4),sprintf('n_b_a_r = %6.3g 1e19 m^-^3',zs.nbar(k)./1e19));
    text(Rmax - ex(3), Zmax - 5 .* ex(4),sprintf('f_g_r = %6.3g ',zs.nbar(k)./negr(k)));
    
    fex = 0.7;
    text(Rmax - ex(3).*fex, Zmin + 5 .* ex(4),sprintf('beta_P = %6.3g ',zs.betap(k)));
    text(Rmax - ex(3).*fex, Zmin + 4 .* ex(4),sprintf('l_i = %6.3g ',zs.li(k)));
    text(Rmax - ex(3).*fex, Zmin + 3 .* ex(4),sprintf('q_a = %6.3g ',profli.qjli(kf,end)));
    text(Rmax - ex(3).*fex, Zmin + 2 .* ex(4),sprintf('q_0 = %6.3g ',profli.qjli(kf,2)));
    text(Rmax - ex(3).*fex, Zmin + 1 .* ex(4),sprintf('q_m_i_n = %6.3g ',min(profli.qjli(kf,:))));
    indqmin = max(find(profli.qjli(kf,:)== min(profli.qjli(kf,:))));
    text(Rmax - ex(3).*fex, Zmin + 0 .* ex(4),sprintf('x_m_i_n = %6.3g ',profli.xli(indqmin)));
    
    % recherche du graphe des profils
    hpf_ = findobj(0,'type','figure','tag','z0profview');
    if ~isempty(hpf_)
        hpf_curs =  findobj(hpf_,'type','uicontrol','tag','curseur');
        if ~isempty(hpf_curs)
            set(hpf_curs,'value',profli.temps(kf));
            zplotprof('curseur',hpf_);
            figure(hf);
        end
    end
    %update
    drawnow
    
    % ajout d'une image
    if ischar(filename)
        F = getframe(hf);
        if exist('VideoWriter')
            writeVideo(vidObj,F);
        else
            mov = addframe(mov,F);
        end
    elseif nargout > 0
        F = getframe(hf);
        mov(kf) =F;
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
    first = 0;
end

% fermeture du fichier
if ischar(filename)
    if exist('VideoWriter')
        close(vidObj);
    else
        mov = close(mov);
    end
else
    filename = mov;
end

set(gca,'zcolor',[0 0 0],'xcolor',[0 0 0],'ycolor',[0 0 0],'color',[1 1 1]);
set(findobj(gcf,'type','text'),'color',[0 0 0]);
set(get(gca,'xlabel'),'color',[0 0 0]);
set(get(gca,'ylabel'),'color',[0 0 0]);
set(get(gca,'title'),'color',[0 0 0]);
set(gcf,'color',[1 1 1]);
map = colormap;
colormap(1-map);

