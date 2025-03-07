liste          = {};
liste{end+1}   = 'ppf/@shot/KK3/CPRF';   
liste{end+1}   = 'ppf/@shot/KK3/TPRF';   
liste{end+1}   = 'ppf/@shot/EFIT/Q';
data         = cgcgetjet(post.z0dinput.shot,liste);
tesh         = data.ppf.KK3.TPRF.data;
ttesh        = data.ppf.KK3.TPRF.t;
rcsh         = data.ppf.KK3.TPRF.x;
Rsh          = data.ppf.KK3.CPRF.data;
%
% Temperature electronique precise
%	
liste                = {};
if rcsh(1) > 9
  eval(['liste{1}      = ''ppf/@shot/KK3/TE',int2str(rcsh(1)),''';'])
else
    eval(['liste{1}      = ''ppf/@shot/KK3/TE0',int2str(rcsh(1)),''';'])
end
for k=2:length(rcsh)
  if rcsh(k) > 9
    eval(['liste{end+1} = ''ppf/@shot/KK3/TE',int2str(rcsh(k)),''';'])
  else
    eval(['liste{end+1} = ''ppf/@shot/KK3/TE0',int2str(rcsh(k)),''';'])
  end
end
for k=1:length(rcsh)
  if rcsh(k) > 9
    eval(['liste{end+1} = ''ppf/@shot/KK3/RC',int2str(rcsh(k)),''';'])
  else
    eval(['liste{end+1} = ''ppf/@shot/KK3/RC0',int2str(rcsh(k)),''';'])
  end
end  
data_te = cgcgetjet(post.z0dinput.shot,liste);

teshf 	       = [];
Rshc  	       = [];
newtemps       = [];
for k=1:length(rcsh)
  if rcsh(k) > 9
    eval(['te = data_te.ppf.KK3.TE',int2str(rcsh(k)),'.data;'])
    eval(['tteshf = data_te.ppf.KK3.TE',int2str(rcsh(k)),'.t;'])
    eval(['re = data_te.ppf.KK3.RC',int2str(rcsh(k)),'.data;'])
    eval(['tRshc = data_te.ppf.KK3.RC',int2str(rcsh(k)),'.t;'])
  else
    eval(['te = data_te.ppf.KK3.TE0',int2str(rcsh(k)),'.data;'])
    eval(['tteshf = data_te.ppf.KK3.TE0',int2str(rcsh(k)),'.t;'])
    eval(['re = data_te.ppf.KK3.RC0',int2str(rcsh(k)),'.data;'])
    eval(['tRshc = data_te.ppf.KK3.RC0',int2str(rcsh(k)),'.t;'])
  end
  if ~isempty(te)
	  if isempty(newtemps)
	    newtemps     = (min(tteshf):1e-4:max(tteshf))';
	  end
	  te	       = interp1(tteshf,te,newtemps,'nearest',NaN);
	  re	       = interp1(tRshc,re,newtemps,'nearest',NaN);
	  tteshf         = newtemps;	      
	  teshf	       = [teshf te];	  
	  Rshc	       = [Rshc re];
 end
end

% retreive plasma center
Rc  = interp1(post.profil0d.temps,post.profil0d.Raxe(:,1),newtemps,'linear','extrap');
dd  = abs(Rshc - Rc * ones(1,size(Rshc,2)));
mask_max = (dd == (min(dd,[],2) * ones(1,size(Rshc,2))));
te_st   = sum(teshf .* mask_max,2) ./ sum(mask_max,2);


% time intervalle
tmin = min(post.zerod.temps);
switch post.z0dinput.shot
case 73221
  tmax = 44.2;
case 73224
  tmax = 44.5;
case 78834
  tmax = 46.6;
case {78842,738842}
  post.z0dinput.shot = 78842;
  tmax = 46.6;
case 79649
  tmax = 44.5;
case 89723
  tmax = 53.4;
otherwise
  tmax = max(post.zerod.temps);
end
indt = find((post.zerod.temps >= tmin) & (post.zerod.temps <= tmax));


h = findobj(0,'type','figure','tag','z0plotstjet');
if isempty(h)
       h=figure('tag','z0plotstjet');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

k    = 2;
subplot(k,1,1);	
plot(post.zerod.temps,medfilt1(post.zerod.te0,3)./1e3,'r',tteshf,te_st./1e3,'b',post.zerod.temps,post.zerod.indice_inv>1,'k');
ylabel('keV')
title(sprintf('%s@%d',post.z0dinput.machine,post.z0dinput.shot));
legend('T_{e,0} METIS','T_{e,0} KK3','METIS S-T indicator');   
set(gca,'xlim',[tmin,tmax]);
subplot(k,1,2);	
plot(post.zerod.temps,post.zerod.q0,'r',post.zerod.temps,post.zerod.qmin,'b',data.ppf.EFIT.Q.t,data.ppf.EFIT.Q.data(:,1),'k', ...
	post.zerod.temps,ones(size(post.zerod.temps)),'g',post.zerod.temps,2 .* ones(size(post.zerod.temps)),'-.g',post.zerod.temps,3/2 .* ones(size(post.zerod.temps)),'-.g');
set(gca,'ylim',[0,ceil(max(post.zerod.q0))])
ylabel('-');
xlabel('time (s)');
legend('q_{0} METIS','q_{min} METIS','q_{0} EFIT');   
set(gca,'xlim',[tmin,tmax]);
edition2
set(findobj(gcf,'type','axes','tag','legend'),'fontsize',16);   