% script de plot de la polarimetrie sur ts
t = data.gene.temps;
if isempty(post)==1

   disp('You need to run the post processing tool')
   return

end
if strcmp(param.from.machine,'TS') | strcmp(param.from.machine,'JET')
nl = post.polar.nl;
nlmes = post.polar.nlmes;
if isempty(nlmes)
  [nnl,tnnl]=tsbase(fix(param.from.shot.num),'snli');
  nlmes = zeros(length(t),5);
  nlmes(:,4) = interp1(tnl,nnl*1e19,t);
end


af = post.polar.af.*180./pi;
afmes = post.polar.afmes.*180./pi;

h1 =findobj(0,'type','figure','tag','polarplotnc1');
if isempty(h1)
     h1=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
               'times','defaultlinelinewidth',3,'tag','polarplotnc1');
else
     figure(h1);
end
clf

cc=get(gca,'colororder');
if isempty(af)
  af = zeros(length(t),5);
end
plot(t,af);
if isempty(afmes)
  afmes = zeros(length(t),5);
end
hold on
set(gca,'colororder',cc);
plot(t,afmes,'o');
st =sprintf('choc %s #%d, t = %4.2g:%4.2g:%4.2g, nbrho = %d',param.from.machine, ...   
             fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
title(st);
grid
xlabel('temps (s)');
ylabel('angle Faraday (degres)');
%set(gca,'xlim',[0 20]);

% correction de l'ecart d'alignement vertical
ind = find(isfinite( prod(af' .* afmes')'));
daf = mean( af(ind,:) -afmes(ind,:));
ind = find(isfinite( prod(nl' .* nlmes')'));
dnl = std(nl(ind,:)-nlmes(ind,:))./mean(nlmes(ind,:));
set(gcf,'ResizeFcn','')

h2 =findobj(0,'type','figure','tag','polarplotnc2');
if isempty(h2)
     h2=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
               'times','defaultlinelinewidth',3,'tag','polarplotnc2');
else
     figure(h2);
end
clf
%figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname','times','defaultlinelinewidth',3);
cc=get(gca,'colororder');
plot(t,af-ones(length(t),1)*daf);
hold on
set(gca,'colororder',cc);
plot(t,afmes,'o');
st =sprintf('choc %s #%d, t = %4.2g:%4.2g:%4.2g, nbrho = %d (correction offset)',param.from.machine, ...              
            fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
title(st);
grid
xlabel('temps (s)');
ylabel('angle Faraday (degres)');
legend(num2str(daf(1),2),num2str(daf(2),2),num2str(daf(3),2),num2str(daf(4),2),num2str(daf(5),2));
%gtext(mat2str(dnl));
%set(gca,'xlim',[0 20]);
set(gcf,'ResizeFcn','')

h3 =findobj(0,'type','figure','tag','polarplotnc3');
if isempty(h3)
     h3=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
               'times','defaultlinelinewidth',3,'tag','polarplotnc3');
else
     figure(h3);
end
clf

% etalonage de la polarimetrie
afcal = NaN * af;
pp=[];
for k =1:size(af,2)
ind = find(isfinite(afmes(:,k))&isfinite(af(:,k)));
%  pp(k,:)      = polyfit(afmes(ind,k),af(ind,k),1);
%  afmesc(:,k)  = polyval(pp(k,:),afmes(:,k));
%  %legstr{k}    = sprintf(' %4.3g \\alpha ^2  %+4.3g \\alpha %+4.3g',pp(k,:));
   pp(k,2) = mean( af(ind,k) -afmes(ind,k));
   afoff   = af(:,k) - pp(k,2);
 	rapp    = afmes(ind,k) ./ afoff(ind)
	indok   = find(abs(afoff(ind)) > 0.2 )
   pp(k,1) = mean(rapp(indok));
   %pp(k,1) = mean(afmes(ind,k) ./ afoff(ind));
   afcal(:,k)   = pp(k,1) .* afoff;
  legstr{k}    = sprintf('%5.3g (\\alpha %+5.3g)',pp(k,:));
  ind2 = find(isfinite(nlmes(:,k))&isfinite(nl(:,k)));
  legstr2{k}   = sprintf('%4.2g',sqrt(mean(((nl(ind2,k)-nlmes(ind2,k))./nlmes(ind2,k)).^2)));
end

subplot(2,1,1)
cc=get(gca,'colororder');
plot(t,afcal);
hold on
set(gca,'colororder',cc);
plot(t,afmes,'o');
st =sprintf('choc %s #%d, t = %4.2g:%4.2g:%4.2g, nbrho = %d (calibration)',param.from.machine, ...   
             fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
title(st);
grid
xlabel('temps (s)');
ylabel('angle Faraday (degres)');

taille = get(gca,'fontsize');
set(gca,'fontsize',12)
legend(gca,legstr,-1);
set(gca,'fontsize',taille)
pos1 =get(gca,'position');

subplot(2,1,2)
set(gca,'colororder',cc);
hold on
plot(t,nl);
hold on
set(gca,'colororder',cc);
plot(t,nlmes,'o');
grid
xlabel('temps (s)');
ylabel('Nl (m^-^3)');
title('Ecart type des Nl')
taille = get(gca,'fontsize');
set(gca,'fontsize',12)
legend(gca,legstr2,-1);
set(gca,'fontsize',taille);
pos2 =get(gca,'position');
pos2(3:4) = pos1(3:4);
set(gca,'position',pos2)
box on
set(gcf,'ResizeFcn','')


end

if strcmp(param.from.machine,'DIIID')
  nl = post.polar.nl;
  nlmes = post.polar.nlmes;
  h1 =findobj(0,'type','figure','tag','interferometry DIIID');
  if isempty(h1)
     h1=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
               'times','defaultlinelinewidth',3,'tag','polarplotnc1');
  else
     figure(h1);
  end
  clf

  cc=get(gca,'colororder');

  if isempty(nl)
    nl = zeros(length(t),4);
  end
  plot(t,nl);
  if isempty(nlmes)
    nlmes = zeros(length(t),4);
  end
  hold on
  set(gca,'colororder',cc);
  plot(t,nlmes/2,'o');
  st =sprintf('shot %s #%d, t = %4.2g:%4.2g:%4.2g, nbrho = %d',param.from.machine, ...   
             fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
  title(st);
  grid
  xlabel('time (s)');
  ylabel('integrated density (m^-^3)');

end
