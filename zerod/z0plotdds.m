% script du  plot des DDS pour TS
if isempty(strmatch(post.z0dinput.machine,'TS','exact'))
  return;
end

zs   = post.zerod;
cons = post.z0dinput.cons;
exp0d  = post.z0dinput.exp0d;
t    = zs.temps;
shot = fix(post.z0dinput.shot);
[gshr,tsh]         = tsbase(shot,'gshr');
[gshte,tsh]      = tsbase(shot,'gshtenv');
if isempty(gshte)
		[gshte,tsh]      = tsbase(shot,'gshte');
end
if ~isempty(gshte)
	indok              = find(all(gshte>=0,1));
	gshr               = mean(gshr(:,indok),1);
	gshte              = gshte(:,indok);
	indok              = max(find(gshr > 2.35 & gshr <= 2.45));
        if isempty(indok)
               d    = abs(gshr - 2.42);
   	       indok              = max(find(d == min(d)));
        end    
	te0                = gshte(:,indok);
end


h = findobj(0,'type','figure','tag','z0plotdds');
if isempty(h)
       h=figure('tag','z0plotdds');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

k    = 2;
subplot(k,1,1);	
plot(t,zs.te0./1e3,'r',t,zs.te0/1e3.*zs.aitb,':r',tsh,te0,'b');
ylabel('Te0 (keV)')
title(sprintf('Zerod : %s@%d/temperature centrale (r -> 0D, b: -> cronos)', ...
             post.z0dinput.machine,post.z0dinput.shot));
   
subplot(k,1,2);	
plot(t,zs.q0,'r',t,zs.qmin,'b', ...
	t,ones(size(t)),'g',t,2 .* ones(size(t)),'-.g',t,3/2 .* ones(size(t)),'-.g');
set(gca,'ylim',[0,min(5,max(zs.q0))])
ylabel('q_0');
xlabel('time (s)');


