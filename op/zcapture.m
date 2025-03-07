% ZCAPTURE  capture une fenetre graphique et en fait une image
%------------------------------------------------------------------
%
% fichier zcapture.m ->  
%
%
% fonction Matlab 5 :
% 
% Cette fonction capture une fenetre de l'interface grahique Cronos
% et cree son image
%
% syntaxe  :
%	zcapture(hf) ;
%
% entrees
%  hf : handle de la fenetre a "capturer"
%
% sorties :
%  				
% fonction ecrite par J-F Artaud, poste 62 15
% version  1.7  du  29/09/2001  
% 
% liste des modifications : 
%
%--------------------------------------------------------------

function zcapture(hf)

%disp('capture d''image inhibee')
%return

if nargin  == 0
   hf = gcf ;
elseif isempty(hf)
   hf =gcf ;
end

name = get(hf,'name') ;
if isempty(name)
	disp('capture impossible, nom absent !')
	return
end
DefFigMode = get(0,'DefaultFigurePaperPositionMode') ;
set(0,'DefaultFigurePaperPositionMode','auto') ;

tag = get(hf,'tag');
if isempty(tag)
	tag = name ;
end
if isempty(name)
   named = tag ;
else
   named = name ;
end

% nom non anbigu
named = sprintf('%s (Capture)',named) ;
set(hf,'name',named) ;

outfile = sprintf('%s/image/%s.png',getappdata(0,'root'),tag) ;
hw      = warndlg('Capture en cours pour la documentation automatique','Patience ...') ;

figure(hf);
disp('on fait l''image')
pause(1)
print(hf,'-dpng',outfile) ;
delete(hw) ;

%htemp = figure;
%zimshow(outfile);
%ButtonName=questdlg('Voulez-vous conserver cette capture ?','Documentation automatique','Oui','Non','Oui');
%if strcmp(ButtonName,'Non')
%	delete(outfile)
%end
%delete(htemp)

set(hf,'name',name) ;
set(0,'DefaultFigurePaperPositionMode',DefFigMode) ;
