% plot JET profiles for publication
function z0plot_jet_nice_profile(post,time)

% extract
zs   = post.zerod;
op0d = post.z0dinput.option;
cons = post.z0dinput.cons;
geo = post.z0dinput.geo;
ts    = zs.temps;

profli = post.profil0d;
xli    = profli.xli;
t      = profli.temps;

% select time 
tmax=max(post.profil0d.temps); 
tmin_round = fix(min(post.profil0d.temps));
if (nargin < 2) || isempty(time)
     % selection of the time slice
    evalin('base','z0plot_reference_post;');
    subplot(3,1,1);
    title('choose 1 time slice for profiles visualization')
    drawnow	
    [tprofi,void] = ginput(1);
else
     tprofi = mean(time);
end

% index in zs and profli
k0d = iround(ts,tprofi);
k1d = iround(t,tprofi); 


% layout
dpos=-20; pos1=256; pos2=247; pos3=364; pos4=425; 
psize1=1; psize10=1; psize2=1; psize3=1; psize4=1; psize5=1; psize6=1;
nfig=0;lls='-'; posfig=0;
% title
tit=['t = ',num2str(tprofi),' s'];

% color order
co = [
         0         0    1.0000
         0    0.5000         0
    1.0000         0         0
         0    0.7500    0.7500
    0.7500         0    0.7500
    0.7500    0.7500         0
    0.2500    0.2500    0.2500
];    

% NE part
%posfig=posfig+dpos;
%nfig=nfig+1;
figuren(nfig); set(gcf,'name','n_e profiles','Position',[pos1+posfig pos2+posfig pos3*psize1 pos4*psize1]);
set(gca,'FontSize',fix(18*psize10),'defaultlinelinewidth',2,'defaultlinelinestyle',lls)
hold on
leg = {};
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
if ~isempty(da_ne.NE.t) && ~isempty(da_ne.DNE.data)
    thrts = da_ne.NE.t * ones(size(da_ne.NE.x));
elseif ~isempty(da_ne.NE.t) && ~isempty(da_ne.NE.data)
    thrts = da_ne.NE.t;
end
if ~isempty(da_ne.RMID.t)
      rhrts = da_ne.RMID.data;
elseif ~isempty(da_ne.NE.t) && ~isempty(da_ne.NE.data)
      rhrts = ones(size(da_ne.NE.t)) * da_ne.NE.x;
else 
      rhrts = [];
end
if ~isempty(da_ne.DNE.t) &&  ~isempty(da_ne.DNE.data)
    khrhts = iround(thrts(:,1),tprofi);
    nehrts=da_ne.NE.data(khrhts,:);
    ne_plus  = da_ne.NE.data(khrhts,:) + da_ne.DNE.data(khrhts,:);
    ne_moins = da_ne.NE.data(khrhts,:) - da_ne.DNE.data(khrhts,:);
    ne_max = 3 * ceil(max(profli.nep(:)./1e19))*1e19;
    indbad = find((ne_moins < 0) | (ne_plus > ne_max));
    ne_plus(indbad) = NaN;
    ne_moins(indbad) = NaN;
    nehrts(indbad) = NaN;
    nehrts_e=[ne_moins(:)';ne_plus(:)';ones(size(da_ne.NE.data(khrhts,:)))*NaN];
    nehrts_e= reshape(nehrts_e(:),length(nehrts_e(:))/length(khrhts),length(khrhts))';
    rhrts_e=[rhrts(khrhts,:);rhrts(khrhts,:);rhrts(khrhts,:)];
    rhrts_e= reshape(rhrts_e(:),length(rhrts_e(:))/length(khrhts),length(khrhts))';
    try
	  plot(rhrts(khrhts,:),nehrts./1e19, 'color',co(1,:),'linestyle','none','marker','o');
	  leg{end+1} = 'n_e HRTS';
    end
elseif ~isempty(da_ne.NE.t) && ~isempty(da_ne.NE.data)
    try
	khrhts = iround(thrts(:,1),tprofi);
	nehrts=da_ne.NE.data(khrhts,:);
	plot(rhrts(khrhts,:),nehrts./1e19,'color',co(1,:),'linestyle','none','marker','o');
	leg{end+1} = 'n_e HRTS';
    end
end

liste        = {};
liste{end+1} = 'ppf/@shot/LIDR/NE';
da           = cgcgetjet(post.z0dinput.shot,liste,'','');
da           = da.ppf.LIDR;
try
    klid = iround(da.NE.t,tprofi);
    plot(da.NE.x,da.NE.data(klid,:)./1e19,'color',co(2,:),'linestyle','none','marker','s');
    leg{end+1} = 'n_e LIDAR';
end
amat = interp1(zs.temps,geo.a,profli.temps,'nearest');
rli   = profli.Raxe + amat * profli.xli;
rmax  = ceil(max(rli(:))*2) / 2;
plot(rli(k1d,:),profli.nep(k1d,:)./1e19,'color',co(3,:));
leg{end+1} = 'n_e METIS';
rli   = profli.Raxe - amat * profli.xli;
rmin  = fix(min(rli(:)));
plot(rli(k1d,:),profli.nep(k1d,:)./1e19,'color',co(3,:));
try
      plot(rhrts_e,nehrts_e./1e19,'color',co(1,:),'linestyle','-','marker','none');         
end
xlabel('R (m)');
ylabel('10^{19} m^{-3}');
legend(leg,'Location','best');
title(tit);
axis([rmin,rmax,0,Inf]);

% Te part
posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','Te profiles','Position',[pos1+posfig pos2+posfig pos3*psize1 pos4*psize1]);
set(gca,'FontSize',fix(18*psize10),'defaultlinelinewidth',2,'defaultlinelinestyle',lls)
hold on

leg  = {};
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
if ~isempty(da_te.TE.t)  &&  ~isempty(da_te.DTE.data)
      thrts = da_te.TE.t * ones(size(da_te.TE.x));
elseif ~isempty(da_te.TE.t) && ~isempty(da_te.TE.data)
    thrts = da_te.TE.t;
end
if ~isempty(da_te.RMID.t)
	rhrts = da_te.RMID.data;
elseif ~isempty(da_te.TE.t)  &&  ~isempty(da_te.DTE.data)
	rhrts = ones(size(da_te.TE.t)) * da_te.TE.x;
else
    rhrts = [];
end
if ~isempty(da_te.DTE.t) &&  ~isempty(da_te.DTE.data)
    khrhts = iround(thrts(:,1),tprofi);
    tehrts=da_te.TE.data(khrhts,:);
    te_plus  = da_te.TE.data(khrhts,:) + da_te.DTE.data(khrhts,:);
    te_moins = da_te.TE.data(khrhts,:) - da_te.DTE.data(khrhts,:);
    te_max = 3 * ceil(max(profli.tep(:)./1e3))*1e3;
    indbad = find((te_moins < 0) | (te_plus > te_max));
    te_plus(indbad) = NaN;
    te_moins(indbad) = NaN;
    tehrts(indbad) = NaN;
    tehrts_e=[te_moins(:)';te_plus(:)';ones(size(da_te.TE.data(khrhts,:)))*NaN];
    tehrts_e= reshape(tehrts_e(:),length(tehrts_e(:))/length(khrhts),length(khrhts))';
    rhrts_e=[rhrts(khrhts,:);rhrts(khrhts,:);rhrts(khrhts,:)];
    rhrts_e= reshape(rhrts_e(:),length(tehrts_e(:))/length(khrhts),length(khrhts))';
    try
        plot(rhrts(khrhts,:),tehrts./1e3, 'color',co(1,:),'linestyle','none','marker','o');
	leg{end+1} = 'T_e HRTS';
  end
else
    try
	khrhts = iround(thrts(:,1),tprofi);
	tehrts=da_te.TE.data(khrhts,:);
	plot(rhrts(khrhts,:),tehrts./1e3,'color',co(1,:),'linestyle','none','marker','o');
	leg{end+1} = 'T_e HRTS';
    end
end
liste        = {};
liste{end+1} = 'ppf/@shot/LIDR/TE';
da           = cgcgetjet(post.z0dinput.shot,liste,'','');
da           = da.ppf.LIDR;
try 
    klid = iround(da.TE.t,tprofi);
    plot(da.TE.x,da.TE.data(klid,:)./1e3,'color',co(2,:),'linestyle','none','marker','s');
    leg{end+1} = 'T_e LIDAR';
end

% ECE
liste        = {};
liste{end+1} = 'ppf/@shot/KK3/TPRF';
liste{end+1} = 'ppf/@shot/KK3/CPRF'; 
da           = cgcgetjet(post.z0dinput.shot,liste,'','');
da           = da.ppf.KK3;
if ~isempty(da.TPRF.data) && ~isempty(da.CPRF.data)
  tesh         = da.TPRF.data;
  ttesh        = da.TPRF.t;
  rcsh         = da.TPRF.x;
  Rsh          = da.CPRF.data;
  kece = iround(ttesh,tprofi);
  try 
    plot(Rsh(kece,:),tesh(kece,:)./1e3,'color',co(4,:),'linestyle','none','marker','d');
    leg{end+1} = 'T_e ECE (KK3)';
  end
end
amat = interp1(zs.temps,geo.a,profli.temps,'nearest');
rli   = profli.Raxe + amat * profli.xli;
rmax  = ceil(max(rli(:))*2) / 2;
plot(rli(k1d,:),profli.tep(k1d,:)./1e3,'color',co(3,:));
leg{end+1} = 'T_e METIS';
rli   = profli.Raxe - amat * profli.xli;
rmin  = fix(min(rli(:)));
plot(rli(k1d,:),profli.tep(k1d,:)./1e3,'color',co(3,:));
try
      plot(rhrts_e,tehrts_e./1e3,'color',co(1,:),'linestyle','-','marker','none');         
end
xlabel('R (m)');
ylabel('keV');
legend(leg,'Location','best');
title(tit);
axis([rmin,rmax,0,Inf]);





% Ti part
posfig=posfig+dpos;
nfig=nfig+1;
figuren(nfig); set(gcf,'name','Ti profiles','Position',[pos1+posfig pos2+posfig pos3*psize1 pos4*psize1]);
set(gca,'FontSize',fix(18*psize10),'defaultlinelinewidth',2,'defaultlinelinestyle',lls)
hold on

leg = {};
liste        = {};
liste{end+1}   = 'ppf/@shot/CXFM/TI';      % profil de Ti
liste{end+1}   = 'ppf/@shot/CXFM/RCOR';      % R of measurment
liste{end+1}   = 'ppf/@shot/CXFM/TILO';    %    eV  minimum  Ti 	
liste{end+1}   = 'ppf/@shot/CXFM/TIHI';    %	eV maximum	 Ti 	
liste{end+1}   = 'ppf/@shot/CXFM/TICR';      % profil de Ti
liste{end+1}   = 'ppf/@shot/CXFM/TIRH';      % profil de Ti
liste = switchuid(liste,post.z0dinput.shot);
da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
da_ti          = da_ti.ppf.CXFM;
if isempty(da_ti.TI.data)
  liste        = {};
  liste{end+1}   = 'ppf/@shot/CXSM/TI';      % profil de Ti
  liste{end+1}   = 'ppf/@shot/CXSM/RCOR';      % R of measurment
  liste{end+1}   = 'ppf/@shot/CXSM/TILO';    %    eV  minimum  Ti 	
  liste{end+1}   = 'ppf/@shot/CXSM/TIHI';    %	eV maximum	 Ti 	
  liste{end+1}   = 'ppf/@shot/CXSM/TICR';      % profil de Ti
  liste{end+1}   = 'ppf/@shot/CXSM/TIRH';      % profil de Ti
  liste = switchuid(liste,post.z0dinput.shot);
  da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
  da_ti          = da_ti.ppf.CXSM;  
end
if isempty(da_ti.TI.data)
  liste        = {};
  liste{end+1}   = 'ppf/@shot/CXGM/TI';      % profil de Ti
  liste{end+1}   = 'ppf/@shot/CXGM/RCOR';      % R of measurment
  liste{end+1}   = 'ppf/@shot/CXGM/TILO';    %    eV  minimum  Ti 	
  liste{end+1}   = 'ppf/@shot/CXGM/TIHI';    %	eV maximum	 Ti 	
  liste{end+1}   = 'ppf/@shot/CXGM/TICR';      % profil de Ti
  liste{end+1}   = 'ppf/@shot/CXGM/TIRH';      % profil de Ti
  liste = switchuid(liste,post.z0dinput.shot);
  da_ti          = cgcgetjet(post.z0dinput.shot,liste,'','');
  da_ti          = da_ti.ppf.CXGM;  
end
if ~isempty(da_ti.TI.t)  && ~isempty(da_ti.TILO.data)
      thrts = da_ti.TI.t * ones(size(da_ti.TI.x));
elseif ~isempty(da_ti.TI.t)  && ~isempty(da_ti.TI.data)
      thrts = da_ti.TI.t;
end
if ~isempty(da_ti.RCOR.t)
	rhrts = da_ti.RCOR.data;
elseif ~isempty(da_ti.TI.t)  && ~isempty(da_ti.TILO.data)
	rhrts = ones(size(da_ti.TI.t)) * da_ti.TI.x;
else
	rhrts = [];
end
if ~isempty(da_ti.TILO.t) && ~isempty(da_ti.TILO.data)
    khrhts = iround(thrts(:,1),tprofi);
    if isfinite(khrhts)
      
	tihrts=da_ti.TI.data(khrhts,:);
	ti_plus  = da_ti.TIHI.data(khrhts,:);
	ti_moins = da_ti.TILO.data(khrhts,:);
	ti_max = 3 * ceil(max(profli.tip(:)./1e3))*1e3;
	indbad = find((ti_moins < 0) | (ti_plus > ti_max));
	ti_plus(indbad) = NaN;
	ti_moins(indbad) = NaN;
	tuhrts(indbad) = NaN;
	tihrts_e=[ti_moins(:)';ti_plus(:)';ones(size(da_ti.TI.data(khrhts,:)))*NaN];
	tihrts_e= reshape(tihrts_e(:),length(tihrts_e(:))/length(khrhts),length(khrhts))';
	rhrts_e=[rhrts(khrhts,:);rhrts(khrhts,:);rhrts(khrhts,:)];
	rhrts_e= reshape(rhrts_e(:),length(tihrts_e(:))/length(khrhts),length(khrhts))';
	try
	    plot(rhrts(khrhts,:),tihrts./1e3, 'color',co(1,:),'linestyle','none','marker','o');
	    leg{end+1} = 'T_i CX';
	end
    end
else
    try
	khrhts = iround(thrts(:,1),tprofi);
	tihrts=da_ti.TI.data(khrhts,:);
	plot(rhrts(khrhts,:),tihrts./1e3,'color',co(1,:),'linestyle','none','marker','o');
	leg{end+1} ='T_i CX';
   end 

end
try
      khrhts = iround(thrts(:,1),tprofi);
      plot(rhrts(khrhts,:),da_ti.TICR.data(khrhts,:)./1e3,'color',co(2,:),'linestyle','none','marker','s');
      leg{end+1} ='T_i CX corrected';
end 
  
amat = interp1(zs.temps,geo.a,profli.temps,'nearest');
rli   = profli.Raxe + amat * profli.xli;
rmax  = ceil(max(rli(:))*2) / 2;
plot(rli(k1d,:),profli.tip(k1d,:)./1e3,'color',co(3,:));
leg{end+1} = 'T_i METIS';
rli   = profli.Raxe - amat * profli.xli;
rmin  = fix(min(rli(:)));
plot(rli(k1d,:),profli.tip(k1d,:)./1e3,'color',co(3,:));
try
      plot(rhrts_e,tihrts_e./1e3,'color',co(1,:),'linestyle','-','marker','none');         
end
xlabel('R (m)');
ylabel('keV');
legend(leg,'Location','best');
title(tit);
axis([rmin,rmax,0,Inf]);



% recover figure handle
function h = figuren(n)

tag = sprintf('z0plot_jet_nice_profile_%d',n);
h = findobj(0,'type','figure','tag',tag);
if isempty(h)
       h=figure('tag',tag,'color',[1 1 1]);
else
       figure(h);
       set(h,'color',[1 1 1]);
end   
clf
setcoloroldcolororder(h);

function setcoloroldcolororder(hfig)

if nargin == 0
  hfig = gcf;
end

co = [
         0         0    1.0000
         0    0.5000         0
    1.0000         0         0
         0    0.7500    0.7500
    0.7500         0    0.7500
    0.7500    0.7500         0
    0.2500    0.2500    0.2500
];    



set(hfig,'defaultAxesColorOrder',co);



function xlabelrho

if verLessThan('matlab','8.6')
  xlabel('r','fontname','symbol');
else
  xlabel('\rho');
end


function liste = switchuid(liste,shot)

switch shot
case 86614
  for k=1:length(liste)
    liste{k} = sprintf('%s?uid=cgiroud+seq=977',liste{k});
  end
end

