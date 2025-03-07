%  INPUTD  courte description  
%------------------------------------------------------------------------------- 
% fichier :  inputd.m  ->  inputd 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   [vallu]=inputd(string,valdef) 
%  
% entrees :  
%  string = 
%  valdef = 
%  
% sorties :  
%   [vallu] = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  3.1  du  20/02/2006  
%  
%@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function [vallu]=inputd(string,valdef)

if isempty(valdef),valdef=nan;end
vallu=nan;
while isnan(vallu)
  if ischar(valdef)
    rep=input([string ' ? [',valdef,'] '],'s');
  else
    rep=input(sprintf([string ' ? [%g] '],valdef));
  end
  if ~isempty(rep),vallu=rep;
  else;vallu=valdef;end
end
