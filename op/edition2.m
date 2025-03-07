function edition2(epaisseur)
% change les fontes de la figure courante en taille 18
% fait l'inversion des couleurs egalement

taille=18;  % taille des caracteres
if nargin == 0
   epaisseur=2;  % epaisseur de trait par defaut
end

set(gcf,'Color',[1 1 1]);
%set(gcf,'InvertHardcopy','off');

h=get(gcf,'children');

% boucle sur les axes eventuellement multiples de la figure
for j=1:length(h)
  try
	nax = h(j);
	set(nax,'FontSize',taille);
	set(nax,'LineWidth',epaisseur);
	child=get(nax,'Children');
	for k=1:length(child)
	if ~strcmp(get(child(k),'type'),'text')
		set(child(k),'LineWidth',epaisseur)
	else	
	%        set(child(k),'Fontsize',18)
	end	
	end
	xlab=get(nax,'Xlabel');
	set(xlab,'FontSize',taille);
	ylab=get(nax,'Ylabel');
	set(ylab,'FontSize',taille);
	titre=get(nax,'Title');
	set(titre,'FontSize',taille);
  end
end
