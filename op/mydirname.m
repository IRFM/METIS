%  MYDIRNAME  courte description  
%------------------------------------------------------------------------------- 
% fichier :  mydirname.m  ->  mydirname 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   dirprint = mydirname(jobpath,charplus) 
%  
% entrees :  
%  jobpath  = 
%  charplus = 
%  
% sorties :  
%   dirprint  = 
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
function dirprint = mydirname(jobpath,charplus)

[dum1,jobname,dum2,dum3] = fileparts(jobpath);
tname   = time_bis;
pidname = getidprocess;

if isempty(charplus)
  dirprint = sprintf('%s_ID.%d_%s',jobname,pidname,tname);
else
  dirprint = sprintf('%s_ID.%d_%s_%s',jobname,pidname,tname,charplus);
end

function y = time_bis()
% TIME          Return time in a string.
%function y = time()
%       S = TIME returns a string containing the time
%           in dd-mm-yy  hh:mm:ss format.

t = round(clock); an=int2str(t(1));
mois = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';...
        'Sep';'Oct';'Nov';'Dec'];

y = ['DATE_',int2str(t(3)),'-',mois(t(2),:),'-',an(3:4),'_TIME_',...
     int2str(t(4)),':',int2str(t(5)),':',int2str(t(6))];



