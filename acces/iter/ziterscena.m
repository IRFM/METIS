%  ZITERSCENA  courte description  
%------------------------------------------------------------------------------- 
% fichier :  ziterscena.m  ->  ziterscena 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   sc = ziterscena(choix) 
%  
% entrees :  
%  choix = 
%  
% sorties :  
%   sc  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  1.9  du  11/03/2002  
%  
%*#auto#   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function sc = ziterscena(choix)

if nargin == 0 
	choix = NaN;
elseif isempty(choix)
	choix = NaN;
end

switch choix
    case 1	
       % donnees scenario standart (modulation des valeurs)
	% temps
	sc.t    = [0     10    200     230     600    700    750    800];
	% ip 
	sc.ip   = [eps     0.05  1       1       1      0.8    0.7      eps];
	% zeff
	sc.zeff = [1.5   0.75  1       1       1      1      1     0.75];
	% pidn
	sc.pidn = [eps     0     0.1     1       1      1      eps        eps];
	% pfci
	sc.pfci = [eps     0     0.1     1       1      1      eps        eps];
	% pfce
	sc.pfce = [eps     0     0.1     1       1      1      eps        eps];
	% fraction de tritium ou hydrogene
	sc.ftri = [eps     0     0       1       1      1      0.3      eps];
	% frction d'helium
	sc.fhe  = [0.1   0.1   0.1     1       1      1      0.5      eps];
	% nbar
	sc.nbar = [0.1     0.2   0.5     1       1      1      0.3      eps];
    case 2	
       % donnees scenario standart (modulation des valeurs)
	% temps
	sc.t    = [0     10    200     230     600    700    750    800];
	% ip 
	sc.ip   = [eps     0.05  1       1       1      0.8    0.7      eps];
	% zeff
	sc.zeff = [1.5   0.75  1       1       1      1      1     0.75];
	% pidn
	sc.pidn = [eps     0     0.1     1       1      1      eps        eps];
	% pfci
	sc.pfci = [eps     0     0.1     1       1      1      eps        eps];
	% pfce
	sc.pfce = [eps     0     0.1     1       1      1      eps        eps];
	% fraction de tritium ou hydrogene
	sc.ftri = [eps     0     0       1       1      1      0.3      eps];
	% frction d'helium
	sc.fhe  = [0.1   0.1   0.1     1       1      1      0.5      eps];
	% nbar
	sc.nbar = [0.1     0.2   0.5     1       1      1      0.3      eps];
	
otherwise
	% donnees scenario standart (modulation des valeurs)
	% temps
	sc.t    = [0     10    200     230     600    700    750    800];
	% ip 
	sc.ip   = [eps     0.05  1       1       1      0.8    0.7      eps];
	% zeff
	sc.zeff = [1.5   0.75  1       1       1      1      1     0.75];
	% pidn
	sc.pidn = [eps     0     0.1     1       1      1      eps        eps];
	% pfci
	sc.pfci = [eps     0     0.1     1       1      1      eps        eps];
	% pfce
	sc.pfce = [eps     0     0.1     1       1      1      eps        eps];
	% fraction de tritium ou hydrogene
	sc.ftri = [eps     0     0       1       1      1      0.3      eps];
	% frction d'helium
	sc.fhe  = [0.1   0.1   0.1     1       1      1      0.5      eps];
	% nbar
	sc.nbar = [0.1     0.2   0.5     1       1      1      0.3      eps];
end
