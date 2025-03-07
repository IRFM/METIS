%  LOAD_TOKAMAK  courte description  
%------------------------------------------------------------------------------- 
% fichier :  load_tokamak.m  ->  load_tokamak 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   [tokamak]=load_tokamak(name) 
%  
% entrees :  
%  name = 
%  
% sorties :  
%   [tokamak] = 
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
function [tokamak]=load_tokamak(name)

xmlfile = ['dina' name '.xml'];
matfile = ['dina' name '.mat'];

reload = 1;matdate = 0;
file = dir(which(matfile));
if(size(file,1)) matdate = datenum(file.date);end
file = dir(which(xmlfile));
xmldate = datenum(file.date);
if(xmldate < matdate) reload=0;end

if(reload)
  disp('Loading from XML file')
  eval(['tokamak=convert(xmltree(''' which(xmlfile) '''));']);
  save(matfile,'tokamak')
else
  disp('Loading from Matlab cache')
  load(matfile)
end

disp('Loading complete')
