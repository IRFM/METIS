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
neutrons             =interp1(RAD,neu2(:,2),rho,'linear')+interp1(RAD,neu3(:,2),rho,'linear');
neutronsth           =interp1(RAD,neu1(:,2),rho,'linear')*1e6;
sortie.neutron.dd    = neutrons*1e6;
temp = load('/usr/drfc/cgc/cgc_data/jet/data/53521/temp53521.mat');
if ~isempty(temp.rnt)
    fnexp = temp.rnt;
    tfnexp = temp.trnt;
end
xfn              = gene.x;
sv               = sigvddn(Ti(:,2))'; 
nd               = n1(:,2)'*1e19;
flux.neutron     = 1/2 * equi.rhomax .* trapz(xfn,equi.vpr .* sv .* nd .* nd,2);

flux.neutrontot  = 1/2 * equi.rhomax .* trapz(xfn,equi.vpr .* sortie.neutron.dd,2);
flux.neutronthe  = 1/2 * equi.rhomax .* trapz(xfn,equi.vpr .* neutronsth,2);
flux.neutronexp  = fnexp(min(find(tfnexp>max(temps))));
