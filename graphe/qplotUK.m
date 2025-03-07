% script de plot du profil d q
chargedata
if exist('temps','var')
   tempsold =temps;
   temps = input(['liste des temps a visualiser', mat2str(tempsold)]);
   if isempty(temps )
         temps =tempsold;
    end
 else
    temps =[];
    while isempty(temps)
          temps = input('liste des temps a visualiser ? ');
    end
 end
ind =iround(t,temps);
if ~isempty(jeux1)
  ind1 = iround(t1,temps);
end
texteleg ={};
for k=1:length(temps)
    texteleg{k}=sprintf('t = %4g (s)',temps(k));
end
    
h = findobj(0,'type','figure','tag','qplot1');
if isempty(h)
       h=figure('tag','qplot1');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1])
cc=get(gca,'colororder');
set(gca,'colororder',cc);
plot(x,data.equi.q(ind,:),'-');
legend(texteleg);
hold on
set(gca,'colororder',cc);
if ~isempty(t1)
      plot(x(1:3:end),jeux1.data.equi.q(ind1,1:3:end),'o');
end
plot([0 1],[1 1],'g:');
plot([0 1],[1.5 1.5],'c:');
plot([0 1],[2 2],'b:');
hold off
xlabel('sqrt(phi/pi/b0) normalise (su)')
ylabel('q')
if ~isempty(t1)
    title(['Shot #',int2str(fix(param.from.shot.num)),', - : heat transport + current diffusion, o : current diffusion']);
else
    title(['Shot #',int2str(fix(param.from.shot.num))]);
end

h = findobj(0,'type','figure','tag','qplot2');
if isempty(h)
       h=figure('tag','qplot2');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1])
cc=get(gca,'colororder');
set(gca,'colororder',cc);
plot(x,data.prof.te(ind,:)/1e3,'-');
legend(texteleg);
hold on
set(gca,'colororder',cc);
if ~isempty(t1)
plot(x(1:3:end),jeux1.data.prof.te(ind1,1:3:end)/1e3,'o');
end
hold off
xlabel('sqrt(phi/pi/b0) normalise (su)')
ylabel('Te (keV)')

if ~isempty(t1)
    title(['Shot #',int2str(fix(param.from.shot.num)),', - : heat transport + current diffusion, o : current diffusion']);
else
    title(['Shot #',int2str(fix(param.from.shot.num))]);
end

h = findobj(0,'type','figure','tag','qplot3');
if isempty(h)
       h=figure('tag','qplot3');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1])
cc=get(gca,'colororder');
plot(x,data.source.hyb.j(ind,:));
johm = data.prof.jmoy-data.source.totale.j;
hold on
plot(x,johm(ind,:),'--');
hold off
if ~isempty(jeux1)
hold on
plot(x,jeux1.data.source.hyb.j(ind1,:),'o');
end
legend(texteleg);
xlabel('sqrt(phi/pi/b0)  normalise (su)')
ylabel('A*m^-^2')
if isempty(jeux1)
title(['Shot #',int2str(fix(param.from.shot.num)),', - J_L_H CRONOS, -- Johm']);
else
title(['Shot #',int2str(fix(param.from.shot.num)),', - J_L_H, CRONOS, -- Johm, o : reference']);
end
h = findobj(0,'type','figure','tag','qplot4');
if isempty(h)
       h=figure('tag','qplot4');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1])
cc=get(gca,'colororder');
set(gca,'colororder',cc);
semilogy(x,data.prof.flux.kean(ind,:)./data.prof.ne(ind,:),'-');
legend(texteleg);
hold on
set(gca,'colororder',cc);
if ~isempty(t1)
    semilogy(x(1:3:end),jeux1.data.prof.flux.kean(ind1,(1:3:end))./jeux1.data.prof.ne(ind1,(1:3:end)),'o');
end
semilogy(x,data.neo.coef.ee(ind,:)./data.prof.ne(ind,:),':');
hold off
xlabel('sqrt(phi/pi/b0) normalise (su)')
ylabel('\Xi_e (m^2s^-^1)')
axis([-inf inf 0 50])
title('Shot #28334, - : heat transport , o : from experimental data, ... : neoclassical ');
if ~isempty(t1)
    title(['Shot #',int2str(fix(param.from.shot.num)),', - : heat transport , o : from experimental data, ... : neoclassical ']);
else
     title(['Shot #',int2str(fix(param.from.shot.num)),', - : heat transport , ... : neoclassical ']);
end


h = findobj(0,'type','figure','tag','qplot5');
if isempty(h)
       h=figure('tag','qplot5');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1])
cc=get(gca,'colororder');
set(gca,'colororder',cc);
plot(x,data.source.totale.el(ind,:)/1e6,'-');
legend(texteleg);
hold on
set(gca,'colororder',cc);
if ~isempty(jeux1)
   plot(x(1:3:end),jeux1.data.source.totale.el(ind1,(1:3:end))/1e6,'o');
end
hold off
xlabel('sqrt(phi/pi/b0) normalise (su)')
ylabel('Source_e_l (MW/m^3)')
title('Shot #28334, - : heat transport , o : from experimental data');
if ~isempty(t1)
    title(['Shot #',int2str(fix(param.from.shot.num)),', - : heat transport , o : from experimental data,']);
else
     title(['Shot #',int2str(fix(param.from.shot.num))]);
end
