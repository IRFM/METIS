%  ZCLEAR0DSEPA  courte description  
%------------------------------------------------------------------------------- 
% fichier :  	zclear0dsepa.m  ->  zclear0dsepa 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   zclear0dsepa 
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
function zclear0dsepa

try
	evalin('base','z0dinput.exp0d = rmfield(z0dinput.exp0d,''Rsepa'');');
	evalin('base','z0dinput.exp0d = rmfield(z0dinput.exp0d,''Zsepa'');');
end