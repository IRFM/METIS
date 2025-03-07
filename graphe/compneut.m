%  NOM DE LA FONCTION  courte description  
%------------------------------------------------------------------------------- 
% fichier : nom_du_fichier -> nom_fonction_principale, nom_fonction_locale 1,... 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function .... 
%  
% entrees :  
%  
%  
% sorties :  
%  
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  3.1  du  18/11/2005  
%  
%@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
indd             =  min(find((param.compo.z == 1) & (param.compo.a == 2)));
sv               = sigvddn(data.prof.ti/1e3);

sv1               = sigvddn(jeux1.data.prof.ti/1e3);

xfn              = param.gene.x;
xxfn             = ones(length(data.gene.temps),1) * xfn;
nd               = squeeze(data.impur.impur(:,:,indd));
nd               = tsplinet(rsa,nd,xxfn);
nd1               = squeeze(jeux1.data.impur.impur(:,:,indd));
nd1               = tsplinet(rsa,nd1,xxfn);
flux.neutron     = 1/2 * data.equi.rhomax .* trapz(xfn,data.equi.vpr .* sv .* nd .* nd,2);
flux.neutron1     = 1/2 * jeux1.data.equi.rhomax .* trapz(xfn,jeux1.data.equi.vpr .* sv .* nd1 .* nd1,2);
h = findobj(0,'type','figure','tag','flux neutrons');
if isempty(h)
       h=figure('tag','flux neutrons');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

[fnexp,tfnexp]=tsbase(fix(param.from.shot.num),'gfluntn');
fnexp         = 10.^fnexp;

plot(data.gene.temps,flux.neutron,jeux1.data.gene.temps,flux.neutron1,tfnexp,smooth(fnexp(:,2),4))
legend('Cronos','Reference','experience')
axis([min(data.gene.temps) max(data.gene.temps) -inf inf])
grid
