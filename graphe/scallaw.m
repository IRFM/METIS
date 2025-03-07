function [o1,o2] = scallaw(loi,ip,rgeo,amin,kappa,nel,bt,meff,plth)
%
% function lois         = scallaw
% function fit          = scallaw(loi,ip,rgeo,amin,kappa,nel,bt,meff,plth)
% function [fit,t,tau]  = scallaw(loi,choc[,racine[,cas]])
% function [c,alpha]    = scallaw(loi)
%
% Arguments d entree
%
%    loi           nom de la loi d'echelle
%
%    ip            courant plasma                    (MA)
%    rgeo          grand rayon plasma                (m)
%    amin          petit rayon plasma                (m)
%    kappa         elongation
%    nel           densite centrale lineique moyenne (e19)
%    bt            champ toroidal                    (T)
%    meff          masse effective                   (uma) 
%    plth          puissance totale perdue           (MW)
%
%    choc          numero du choc
%    racine        repertoire ou chercher les donnees par defaut dans la base TS
%    cas           sous-repertoire
%
% Nota Bene
%
%    Coefficient des lois de scaling du temps de confinement thermique
% exprime en secondes de la forme fit = c * (IP ^ alpha(1)) * ... * (PLTH ^alpha(8)) :
%
% c    coefficient multiplicatif
%
% alpha(1) exposant IP   (MA)
% alpha(2) exposant RGEO (m)
% alpha(3) exposant AMIN (m)
% alpha(4) exposant KAPPA
% alpha(5) exposant NEL  (e19)
% alpha(6) exposant BT   (Pa)
% alpha(7) exposant MEFF
% alpha(8) exposant PLTH (MW)
%
% Les lois disponibles sont :
%
%      'ts'         Loi de scaling TS
%      'iterth96'   Loi de scaling ITER THERMIQUE 96
%      'iter0'      Loi de scaling ITER modifiee
%      'iterth93'   Loi de scaling ITER93-H ELM-free
%      'iterth93y'  Loi de scaling ITER93-H ELMy = 0.85 * ITER93-H ELM-free
%      'textorRI'   Loi du mode RI de TEXTOR
%      'tsRI'       Loi du mode RI-like de TS (chauffage minoritaire)
%      'tsEH'       Loi de scaling TS en mode ameliore chauffage aux electrons
%      'textorRIa'  Loi du mode RI de TEXTOR corrigee du AMIN ( = 0.46 m) et 
%                   du courant IP ( = 0.4) et de l'effet de masse 
%
%  'iter0' est une correction de la loi 'iterth96' dans laquelle on a supprime
%  la dependance en MEFF. En effet la loi 'iterth96' colle bien avec la loi 'ts'
%  en Deuterium mais mal en Helium
 
lois = str2mat('iter0','iterth96','ts','iterth93','iterth93y','textorRI');

if nargin == 0

   o1 = lois;
   return;

end

if strcmp(loi,'iterth96')

   c     = 0.023;
   alpha = [0.96 1.89 -0.06 0.64 0.4 0.03 0.2 -0.73];

elseif strcmp(loi,'ts')

   c     = 0.0227;
   alpha = [0.9821 1.8366 0.0 0.0 0.4290 0.2032 0.0 -0.7475];

elseif strcmp(loi,'iter0')

   c     = 0.023 * (2 ^ 0.2);
   alpha = [0.96 1.89 -0.06 0.64 0.4 0.03 0.0 -0.73];

elseif strcmp(loi,'iterth93')

   c     = 0.036;
   alpha = [1.06 1.90 -0.11 0.66 0.17 0.32 0.41 -0.67];

elseif strcmp(loi,'iterth93y')

   c     = 0.85 * 0.036;
   alpha = [1.06 1.90 -0.11 0.66 0.17 0.32 0.41 -0.67];

elseif strcmp(loi,'textorRI')

   c     = 0.0113;
   alpha = [0.06 1.90  1.89 0.66 1.17 0.32 0.41 -0.67];

elseif strcmp(loi,'tsRI')

   c     = 0.1633;
%  alpha = [0.6538 1.90 -0.11 0.66 0.0286 0.32 0.0461 -0.445];
   alpha = [0.6538 0     0    0    0.033 0.24 0      -0.45];

elseif strcmp(loi,'tsEH')

%
% Le coefficient de MEFF a ete impose a zero
% La loi est construite a partir de choc en He seulement
%

   c     = 0.0161 * (4 ^ 0.41);
   alpha = [0.7433 1.90 -0.11 0.66 0.1919 0.3221 0.00 -0.6343];

elseif strcmp(loi,'textorRIa')

   c     = 0.0113 * (0.46 ^ 2) / 0.4 * (4 ^ 0.41);
   alpha = [1.06 1.90 -0.11 0.66 1.17 0.32 0.00 -0.67];

else

   disp('Loi existantes :');
   disp(lois);
   return;
   
end

if nargin >= 2 & nargin <= 4

   choc = ip;
   
   if nargin < 4
   
      cas = '';
      
      if nargin < 3
      
         racine = cgcdefault;

      else
      
         racine = rgeo;
                  
      end
      
   else
   
      cas    = amin;
      racine = rgeo;
         
   end

   [file,etat,mess] = cgcdata(choc,'bile','r',racine,cas);
   
   if etat ~= 1
   
      disp(['scallaw - ' mess]); 
      return;
      
   end

   [ip,rgeo,amin,kappa,nel,bt,mmain,plth] = cgcreadfile(file,...
   'ip','rmaj','amin','elong','nbar','btor','mmain','ploss');
   meff = mmain(1);
           
end         

if nargin > 1

      
%      o1     = c * (ip    .^ alpha(1))...
%                .* (rgeo  .^ alpha(2))...
%                .* (amin  .^ alpha(3))...
%                .* (kappa .^ alpha(4))...
%                .* (nel   .^ alpha(5))...
%                .* (bt    .^ alpha(6))...
%                .* (meff  .^ alpha(7))...
%                .* (plth  .^ alpha(8));

   if length(meff) == 1
   
      meff = ones(size(ip)) * meff;
   
   end

   mask           = alpha ~= 0;
   P              = [ip rgeo amin kappa nel bt meff plth];
   o1             = c * exp(log(P(:,mask)) * alpha(mask)');

   if ~isempty(find(P < 0))
   
      disp('scallaw - some parameter have negative value');
      o1 = real(o1);
   
   end
   
else

   o1 = c;
   o2 = alpha;
   
end
