% constante physique (phys)
phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)


% script pour le plot des donnees HRTS de JEt dans les simulation avec METIS
h = findobj(0,'type','figure','tag','z0plot_jet_edge1');
if isempty(h)
       h=figure('tag','z0plot_jet_edge1');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
colormap('hot')

% NE
subplot(2,2,1)
  
t   = post.profil0d.temps;
rli = post.profil0d.Raxe + interp1(post.z0dinput.cons.temps,post.z0dinput.geo.a,t,'linear','extrap') * post.profil0d.xli;
zplotprof(gca,t,rli,post.profil0d.nep./1e19,'color','r');
%rli   = data.equi.raxe - data.equi.a;
%zplotprof(gca,t,rli,data.prof.ne./1e19,'color','r');

liste        = {};
liste{end+1} = 'ppf/@shot/HRTS/NE';
liste{end+1} = 'ppf/@shot/HRTS/DNE';
liste{end+1} = 'ppf/@shot/HRTS/RMID';     
da_ne          = cgcgetjet(post.z0dinput.shot,liste,'','');
da_ne          = da_ne.ppf.HRTS;
if isempty(da_ne.NE.t) || ischar(da_ne.NE.data) 
	liste        = {};
	liste{end+1} = 'ppf/@shot/HRTS/NE?uid=chain1';
	liste{end+1} = 'ppf/@shot/HRTS/DNE?uid=chain1';
	liste{end+1} = 'ppf/@shot/HRTS/RMID?uid=chain1';     
	da_ne           = cgcgetjet(post.z0dinput.shot,liste,'','');
	da_ne           = da_ne.ppf.HRTS;
end
if isempty(da_ne.NE.t) || ischar(da_ne.NE.data) 
	liste        = {};
	liste{end+1} = 'ppf/@shot/LIDR/NE';
	liste{end+1} = 'ppf/@shot/LIDR/DNE';
	liste{end+1} = 'ppf/@shot/LIDR/RMIDN';   
	da_ne           = cgcgetjet(post.z0dinput.shot,liste,'','');
	da_ne           = da_ne.ppf.LIDR;
    da_ne.RMID      = da_ne.RMIDN;
end
thrts = da_ne.NE.t(:) * ones(1,length(da_ne.NE.x));
if ~isempty(da_ne.RMID.t)
	rhrts = da_ne.RMID.data;
else
	rhrts = size(ones(da_ne.NE.t)) * da_ne.NE.x;
end
if ~isempty(da_ne.DNE.t)
   nehrts=da_ne.NE.data;
   nehrts_e=[da_ne.NE.data(:)' - da_ne.DNE.data(:)';da_ne.NE.data(:)' + da_ne.DNE.data(:)';ones(size(da_ne.NE.data(:)'))*NaN];
   nehrts_e= reshape(nehrts_e(:),length(nehrts_e(:))/length(da_ne.NE.t),length(da_ne.NE.t))';
   rhrts_e=[rhrts(:)';rhrts(:)';rhrts(:)'];
   rhrts_e= reshape(rhrts_e(:),length(rhrts_e(:))/length(da_ne.NE.t),length(da_ne.NE.t))';
   thrts_e=[thrts(:)';thrts(:)';thrts(:)'];
   thrts_e= reshape(thrts_e(:),length(thrts_e(:))/length(da_ne.NE.t),length(da_ne.NE.t))';
   try
	zplotprof(gca,thrts,rhrts,nehrts./1e19, ...
	'color','b','linestyle','none','marker','o');
	zplotprof(gca,thrts_e,rhrts_e,nehrts_e./1e19, ...
	'color','b','linestyle','-','marker','none');
   end
else
    nehrts=da_ne.NE.data
    try
	zplotprof(gca,thrts,rhrts,nehrts./1e19, ...
	'color','b','linestyle','none','marker','o');
    end
end
leg = {'METIS','HRTS'} ;
legend(leg);     
xlabel('R_{ext} (m)')
ylabel('n_e (10^{19} m^{3})');
set(gca,'ylim',[0 Inf]);

% Te
subplot(2,2,2)
%rli   = data.equi.raxe + data.equi.a;
%t = data.gene.temps;
zplotprof(gca,t,rli,post.profil0d.tep./1e3,'color','r');
%rli   = data.equi.raxe - data.equi.a;
%zplotprof(gca,t,rli,data.prof.te./1e3,'color','r');

liste        = {};
liste{end+1} = 'ppf/@shot/HRTS/TE';
liste{end+1} = 'ppf/@shot/HRTS/DTE';
liste{end+1} = 'ppf/@shot/HRTS/RMID';     
da_te          = cgcgetjet(post.z0dinput.shot,liste,'','');
da_te          = da_te.ppf.HRTS;
if isempty(da_te.TE.t) || ischar(da_te.TE.data) 
	liste        = {};
	liste{end+1} = 'ppf/@shot/HRTS/TE?uid=chain1';
	liste{end+1} = 'ppf/@shot/HRTS/DTE?uid=chain1';
	liste{end+1} = 'ppf/@shot/HRTS/RMID?uid=chain1';     
	da_te           = cgcgetjet(post.z0dinput.shot,liste,'','');
	da_te           = da_te.ppf.HRTS;
end
if isempty(da_te.TE.t) || ischar(da_te.TE.data) 
	liste        = {};
	liste{end+1} = 'ppf/@shot/LIDR/TE';
	liste{end+1} = 'ppf/@shot/LIDR/DTE';
	liste{end+1} = 'ppf/@shot/LIDR/RMID';   
	da_te           = cgcgetjet(post.z0dinput.shot,liste,'','');
	da_te           = da_te.ppf.LIDR;
end
thrts = da_te.TE.t(:) * ones(1,length(da_te.TE.x));
if ~isempty(da_te.RMID.t)
	rhrts = da_te.RMID.data;
else
	rhrts = size(ones(da_te.TE.t)) * da_te.TE.x;
end
if ~isempty(da_te.DTE.t)
   tehrts=da_te.TE.data;
   tehrts_e=[da_te.TE.data(:)' - da_te.DTE.data(:)';da_te.TE.data(:)' + da_te.DTE.data(:)';ones(size(da_te.TE.data(:)'))*NaN];
   tehrts_e= reshape(tehrts_e(:),length(tehrts_e(:))/length(da_te.TE.t),length(da_te.TE.t))';
   rhrts_e=[rhrts(:)';rhrts(:)';rhrts(:)'];
   rhrts_e= reshape(rhrts_e(:),length(rhrts_e(:))/length(da_te.TE.t),length(da_te.TE.t))';
   thrts_e=[thrts(:)';thrts(:)';thrts(:)'];
   thrts_e= reshape(thrts_e(:),length(thrts_e(:))/length(da_te.TE.t),length(da_te.TE.t))';
   try
	zplotprof(gca,thrts,rhrts,tehrts./1e3, ...
	'color','b','linestyle','none','marker','o');
	zplotprof(gca,thrts_e,rhrts_e,tehrts_e./1e3, ...
	'color','b','linestyle','-','marker','none');
   end
else
    nehrts=da_te.TE.data
    try
	zplotprof(gca,thrts,rhrts,tehrts./1e3, ...
	'color','b','linestyle','none','marker','o');
    end
end
leg = {'METIS','HRTS'} ;
legend(leg);     
xlabel('R_{ext} (m)')
ylabel('T_e (keV)');
set(gca,'ylim',[0 Inf]);



% Pe
subplot(2,2,3)
%rli   = data.equi.raxe + data.equi.a;
%t = data.gene.temps;
zplotprof(gca,t,rli,(post.profil0d.tep .* post.profil0d.nep) .* phys.e ./1e3,'color','r');
leg = {'METIS','HRTS'} ;
%rli   = data.equi.raxe - data.equi.a;
%zplotprof(gca,t,rli,data.prof.pe./1e3,'color','r');
try 
	zplotprof(gca,thrts,rhrts,da_te.TE.data .* da_ne.NE.data .* phys.e ./ 1e3, ...
	'color','b','linestyle','none','marker','o');
end
if ~isempty(da_te.DTE.t) 
   pehi = (da_te.TE.data  .* da_ne.NE.data + da_te.DTE.data .* da_ne.NE.data + da_te.TE.data .* da_ne.DNE.data) .* phys.e ./ 1e3;
   pelo = (da_te.TE.data  .* da_ne.NE.data - da_te.DTE.data .* da_ne.NE.data - da_te.TE.data .* da_ne.DNE.data) .* phys.e ./ 1e3;
   pehrts_e=[pelo(:)';pehi(:)';ones(size(pelo(:)'))*NaN];
   pehrts_e= reshape(pehrts_e(:),length(pehrts_e(:))/length(da_te.TE.t),length(da_te.TE.t))';
   try
	zplotprof(gca,thrts_e,rhrts_e,pehrts_e, ...
	'color','b','linestyle','-','marker','none');
   end
end
legend(leg);     
xlabel('R_{ext} (m)')
ylabel('P_e (Pa)');
set(gca,'ylim',[0 Inf]);

% Ti
liste        = {};
liste{end+1}   = 'ppf/@shot/CXFM/TIMX';    % Ti0
liste{end+1}   = 'ppf/@shot/CXFM/TI';      % profil de Ti
liste{end+1}   = 'ppf/@shot/CXFM/TICR';      % profil de Ti
liste{end+1}   = 'ppf/@shot/CXFM/TIRH';      % profil de Ti
liste{end+1}   = 'ppf/@shot/CXFM/RCOR';      % R of measurment
liste{end+1}   = 'ppf/@shot/CXFM/TILO';    %    eV  minimum  Ti 	
liste{end+1}   = 'ppf/@shot/CXFM/TIHI';    %	eV maximum	 Ti 	

da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
da_ti          = da_ti.ppf.CXFM;



subplot(2,2,4);

%rext  = data.equi.raxe + data.equi.a;
%t = data.gene.temps;
zplotprof(gca,t,rli,post.profil0d.tip ./ 1e3,'color','r');
leg = {'METIS'} ;
thrts = da_ti.TI.t * ones(size(da_ti.TI.x));
if ~isempty(da_ti.RCOR.t)
	rhrts = da_ti.RCOR.data;
elseif ~isempty(da_ti.TI.t)
	rhrts = size(ones(da_ti.TI.t)) * da_ti.TI.x;
else
	rhrts = [];
end
try
	zplotprof(gca,thrts,rhrts,da_ti.TI.data ./1e3,'color','b','marker','o','linestyle','none');
	leg{end+1} ='CXSM';
end 
try
	zplotprof(gca,thrts,rhrts,da_ti.TICR.data ./1e3,'color','c','marker','+','linestyle','none');
	leg{end+1} ='CXSM corrected';
end 

% plot error bar 
if ~isempty(da_ti.TILO.t)
	tihrts=da_ti.TI.data;
	tihrts_e=[da_ti.TILO.data(:)';da_ti.TIHI.data(:)' ;ones(size(da_ti.TI.data(:)'))*NaN];
	tihrts_e= reshape(tihrts_e(:),length(tihrts_e(:))/length(da_ti.TI.t),length(da_ti.TI.t))';
	rhrts_e=[rhrts(:)';rhrts(:)';rhrts(:)'];
	rhrts_e= reshape(rhrts_e(:),length(rhrts_e(:))/length(da_ti.TI.t),length(da_ti.TI.t))';
	thrts_e=[thrts(:)';thrts(:)';thrts(:)'];
	thrts_e= reshape(thrts_e(:),length(thrts_e(:))/length(da_ti.TI.t),length(da_ti.TI.t))';
	try
		zplotprof(gca,thrts,rhrts,tihrts./1e3, ...
		'color','b','linestyle','none','marker','o');
		zplotprof(gca,thrts_e,rhrts_e,tihrts_e./1e3, ...
		'color','b','linestyle','-','marker','none');
	end
end
legend(leg);     
xlabel('R_{ext} (m)')
ylabel('T_i (keV)');
set(gca,'ylim',[0 Inf]);


% target data : Langmuir probes outer in iner bord
liste        = {};
liste{end+1}   = 'ppf/@shot/Y4PO/NES?uid=KY4D';    % Ti0
liste{end+1}   = 'ppf/@shot/Y4PO/TES?uid=KY4D';      % profil de Ti
liste{end+1}   = 'ppf/@shot/Y4PI/NES?uid=KY4D';    % Ti0
liste{end+1}   = 'ppf/@shot/Y4PI/TES?uid=KY4D';      % profil de Ti

target         = cgcgetjet(post.z0dinput.shot,liste,'','');

if ~isempty(target.ppf.Y4PO.TES.data)
    
    h = findobj(0,'type','figure','tag','z0plot_jet_edge2');
    if isempty(h)
        h=figure('tag','z0plot_jet_edge2');
    else
        figure(h);
    end
    clf
    set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
        'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
    colormap('hot')
    
    
    subplot(2,1,1)
    plot(target.ppf.Y4PO.TES.t,sgolayfilt(target.ppf.Y4PO.TES.data,1,7));
    hold on
    plot(post.zerod.temps,post.zerod.telim,'k');
    ylabel('T_e (eV)');
    title('color = Langmuir, black = METIS (outer target)');
    
    subplot(2,1,2)
    plot(target.ppf.Y4PO.NES.t,sgolayfilt(target.ppf.Y4PO.NES.data,1,7) ./ 1e19);
    hold on
    plot(post.zerod.temps,post.zerod.nelim./ 1e19,'k');
    ylabel('N_e (1e19 m^{-3})');
    xlabel('time (s)');
    set(gca,'ylim',[0,max(target.ppf.Y4PO.NES.data(:))/1e19]);
    
end


if ~isempty(target.ppf.Y4PI.TES.data)
    
    h = findobj(0,'type','figure','tag','z0plot_jet_edge3');
    if isempty(h)
        h=figure('tag','z0plot_jet_edge3');
    else
        figure(h);
    end
    clf
    set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
        'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
    colormap('hot')
    
    
    subplot(2,1,1)
    plot(target.ppf.Y4PI.TES.t,sgolayfilt(target.ppf.Y4PI.TES.data,1,7));
    hold on
    plot(post.zerod.temps,post.zerod.telim,'k');
    ylabel('T_e (eV)');
    title('color = Langmuir, black = METIS (Inner target)');
    
    subplot(2,1,2)
    plot(target.ppf.Y4PI.NES.t,sgolayfilt(target.ppf.Y4PI.NES.data,1,7) ./ 1e19);
    hold on
    plot(post.zerod.temps,post.zerod.nelim./ 1e19,'k');
    ylabel('N_e (1e19 m^{-3})');
    xlabel('time (s)');
    set(gca,'ylim',[0,max(target.ppf.Y4PI.NES.data(:))/1e19]);
    
end
