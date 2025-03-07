% cree les figures pour la comparaison des plots et sauve en png
function h = zplotstruct(s1,s2,nomroot,titre)

if nargin < 3
    nomroot ='';
end
if nargin < 4
	titre = 'noname';
end
h = [];

if isstruct(s1)
	noms = sort(fieldnames(s1));
	nbn  = length(noms);
	
	h = figure;
	clf
	set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
			'defaultlinelinewidth',1,'color',[1 1 1])
	if nargin > 3
		set(h,'name',titre)
	end
	n    = ceil(sqrt(nbn));
	m    = ceil(nbn./n);
	for k = 1:length(noms)
		subplot(n,m,k)
		var1 = getfield(s1,noms{k});
		var1 = var1(:)';
		if isnumeric(var1) & isfield(s2,noms{k}) & (length(size(var1)) == 2)
			var2 = getfield(s2,noms{k});
			var2 = var2(:)';
			if isnumeric(var2) & (length(size(var2)) == 2) & (~isempty(var2)) 
				x1 = linspace(0,1,size(var1,2));
				x2 = linspace(0,1,size(var2,2));
				figure(h);
				if iscomplex(var1)
					plot(x1,real(var1),'ob',x2,real(var2),'r',x1,imag(var1),'oc',x2,imag(var2),'m');
				else
					plot(x1,real(var1),'ob',x2,real(var2),'r');	
				end
				title(noms{k});
			end
			
		set(gca,'ButtonDownFcn','zdataplot(''extrait'');');
		elseif isstruct(var1) & isfield(s2,noms{k})
			var2 = getfield(s2,noms{k});
			zplotstruct(var1,var2,sprintf('%s_%s',nomroot,noms{k}),titre)
			figure(h);
		end
	end
	%eval(sprintf('print(gcf,''-dpng'',''%s'')',nomroot));
elseif isnumeric(s1)  & (length(size(s1)) == 2) & isnumeric(s2)  & (length(size(s2)) == 2)
	h = findobj(0,'type','figure','tag','zplotstruct');
	if isempty(h)
		h=figure('tag','zplotstruct');
	else
		figure(h);
	end   
	clf
	set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
		'defaultlinelinewidth',1,'color',[1 1 1])
	x1 = linspace(0,1,size(s1,2));
	x2 = linspace(0,1,size(s2,2));
	plot(x1,s1,'ob',x2,s2,'r');
	if nargin > 3
		title(titre)
	end
	%eval(sprintf('print -dpng  %s',nomroot));
end