%  ZUICLOSE  remplace delete
%------------------------------------------------------------------------------- 
% fichier : zuicloseone.m  
% 
% 
% fonction Matlab 5 : 
% 
% ferme une fenetre
%  
% syntaxe :  
%   zuiclose
%  
% entrees :  
%  
% sorties :  
%  
% fonction �rite par J-F Artaud , poste 46-78
% version  3.0  du  28/07/2005
%  
%  20/05/2003 -> protection des popup
%  20/05/2003 -> protection objet si ce n'est pas des figures
%  28/07/2005 -> retour a close simple pour essayer sur saturne
%  
% -------------------------------------------------------------------------------
%
function zuicloseone(h)

set(h,'visible','off');
drawnow
delete(h);
drawnow
return


if isempty(h)
   return
end

cp7 = computer;
if strncmp('7',cp7,1)
	close(h);
	return
end

   
if ~checkfigs(h)
    delete(h);
    return
end

hp = findobj(h,'type','uicontrol','style','popupmenu');
if ~isempty(h)
   set(hp,'value',1);
end

if strcmp(computer,'GLNX86')
	if strcmp(get(h,'type'),'figure')
				figure(h);
	end
	set(h,'visible','on')
	pause(0.1);
	drawnow
end 
delete(h) ;
   
%------------------------------------------------
function status = checkfigs(h)
status = 1;
if ~ishandle(h) | ~strcmp(get(h,'type'),'figure')
    status = 0;
    return
end
