% ZUIPLOTIN  substitue des axes a un objet
%-----------------------------------------
% fichier zuiplotin.m ->  zuiplotin
%
% fonction Matlab 5 :
%	Cette fonction  substitue des axes a un objet. 
%
% syntaxe  :
%   ha = zuiplotin(handle,varargin);
	
% entrees :
%	handle    = handle de l'objet
%	varargin  = autres  porprietes des axes
%
% sortie : 
%	ha    = nouvel handle (des axes)  
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 09/02/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function ha=zuiplotin(hui,varargin)

% recuperation des infos
pos    = get(hui,'position');
tag    = get(hui,'tag');
parent = get(hui,'parent');
cb     = get(hui,'callback');
st     = get(hui,'string');
units  = get(hui,'units');

%creation des axes (pour la mesure de l'encombrement
etatmem = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
ha=axes('position',pos,'units',units,'color',[1 1 1],'visible','off', ...
        'tag',tag,'buttondownfcn',cb,'hittest','on'); 
if nargin >1
	set(ha,varargin{:});
end
        
% modification de la liste des handles        
zh = getappdata(parent,'zhandle');
zh = setfield(zh,tag,ha);
setappdata(parent,'zhandle',zh);
% titre et reservation de place
if isempty(st)| all(st ==sprintf(' '))
   st='Ici le titre';
end
   
title(st);
xlabel('Mesure');
ylabel('Mesure');

% calcul de la nouvelle position
drawnow
set(ha,'units','normalized');
drawnow
pos     = get(ha,'position');
extg    = get(get(ha,'ylabel'),'extent');
xgauche = extg(1);
exth    = get(get(ha,'title'),'extent');
yhaut   = exth(2) + exth(4);
extb    = get(get(ha,'xlabel'),'extent');
ybas    = extb(2);
xdroit  = 1.03;

dx = xdroit - xgauche;
dy = (yhaut  - ybas) .*1.03;

newpos(3) = pos(3) ./ dx;
newpos(4) = pos(4) ./ dy;
newpos(1) = pos(1) - newpos(3) .* xgauche;
newpos(2) = pos(2) - newpos(4) .* ybas;


% dessin a la nouvelle position
axes(ha);
xlabel('');
ylabel('');
set(ha,'visible','on','position',newpos); 
zuicloseone(hui);
drawnow
set(ha,'units',units);
set(0,'ShowHiddenHandles',etatmem);
