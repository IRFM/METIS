% script de plot de la polarimetrie sur ts

t =data.gene.temps;
nl = post.polar.nl;
nlmes = post.polar.nlmes;
af = post.polar.af.*180./pi;
afmes = post.polar.afmes.*180./pi;

% correction de l'ecart d'alignement vertical
ind = find(isfinite( prod(af' .* afmes')'));
daf = mean( af(ind,:) -afmes(ind,:));
ind = find(isfinite( prod(nl' .* nlmes')'));
dnl = std(nl(ind,:)-nlmes(ind,:))./mean(nlmes(ind,:));
h    = figure;
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1])


cc=get(gca,'colororder');
plot(t,af-ones(size(t))*daf);
hold on
set(gca,'colororder',cc);
plot(t,afmes,'o');


grid
xlabel('times (s)');
ylabel('Faradays angle (degrees)');

%  legend(['error ='num2str(daf(1),2)],['error ='num2str(daf(2),2)],['error ='num2str(daf(3),2)],...
%         ['error ='num2str(daf(4),2)],['error ='num2str(daf(5),2)],['error ='num2str(daf(6),2)],...
%         ['error ='num2str(daf(7),2)],['error ='num2str(daf(8),2)],-1);
