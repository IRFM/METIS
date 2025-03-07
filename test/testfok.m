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
zineb_compile fokker2tempsnew.f
load testfokker_deneb
      [p1,p2,p3,p4,pel,plo,J,Wp,Wt,WsWth,ds1,ds2,ds3,ds4,ds5,ds6,ds7,ds8,ds9,ds10,ds11,dist,xv,ev,neu1,neu2,neu3,neu4,neu5,neu6,...
      consd1,consd2,encon1,encon2,engkt1,engkt2,consp1,consp2,...
      s2d_1,s2d_2,s2d_3,s2d_4,s2d_5,s2d_6,s2d_7,s2d_8,s2d_9,s2d_10,s2t,...
      RMAT1,RMAT2,RMAT3,RMAT4,RMAT5,RMAT6,RMAT7,RMAT8,RMAT9,RMAT10,RMAT11,fvth,fv,mvpart,mvtot,sourmov,parmom,torque,ndt] = ...
      fokker2tempsnew(geom,beam,XV,plasma,NTVEC,input,pinj1,pinj2,...
                    densite1,densite2,densite3,densite4,densite0,tempion,tempelec,...
                    dep1,dep2,dep3,dep4,dep5,dep6,RAD,neutre,ang,1,vpara,RRz,vol,opt);
