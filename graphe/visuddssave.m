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

time = input('Time ? (s) ');

[tex,t,r] = tsbase(fix(param.from.shot.num),'gshte');
[rte,t,r] = tsbase(fix(param.from.shot.num),'gshr');

[aa,itimeSH] = min(abs(t-time));
icentre = iround(rte(itimeSH,:),2.4); % corde la plus proche du centre

figure
icentf = min(icentre+9,size(tex,2));
nfin   = icentf-icentre;
plot(t,tex(:,[icentre:icentf]))
texfv=[];
for k=icentre:(icentre+nfin)
    if k ~= (icentre+nfin)
    textv = ['num2str(rte(itimeSH,',int2str(k),'),3),'];
    else
    textv = ['num2str(rte(itimeSH,',int2str(k),'),3)'];       
    end
    texfv = cat(2,texfv,textv);
end
legend(num2str(rte(itimeSH,icentre)),num2str(rte(itimeSH,icentre+1)),num2str(rte(itimeSH,icentre+2)),num2str(rte(itimeSH,icentre+3)),num2str(rte(itimeSH,icentre+4)),num2str(rte(itimeSH,icentre+5)),num2str(rte(itimeSH,icentre+6)),num2str(rte(itimeSH,icentre+7)),num2str(rte(itimeSH,icentre+8)),num2str(rte(itimeSH,icentre+9)))

Rinv = input('Enter Rinv as deduced from the figure (m) ');

icronos = iround(data.gene.temps,time);
rhoinv = interp1(squeeze(double(data.equi.R(icronos,:,1))),squeeze(double(data.equi.rhoRZ(icronos,:)))./data.equi.rhomax(icronos),Rinv);

figure
plot(param.gene.x,data.prof.q(icronos,:),[0 1],[1 1],'m--',[rhoinv rhoinv],[0 max(data.prof.q(icronos,:))],'k')
legend('cronos q-profile','q = 1','inv. radius')
xlabel('\rho')
ylabel('q')
