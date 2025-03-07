% ZUIFAITFUNACTION callback de zuicreefunform
%-------------------------------------------------------------------------------------
% fichier zuifaitfunaction.m ->  zuifaitfunaction
%
%
% fonction Matlab 5 :
% 
% Cette fonction gere les callback de zuicreefunform . 
%
% syntaxe  :
%  
%  zuifaitfunaction(action)
%
% entrees :
% 
%  action      = action associee au callback
%
% sortie : aucune
%
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.0, du 27/11/2002.
% 
% 
% liste des modifications : 
%
% * 27/11/2002 -> ajout du code_retour 
%
%--------------------------------------------------------------
%
function zuifaitfunaction(action,hfig)

if nargin ==0
	return;
elseif isempty(action)
	return;
end

% recupere le handle de la fenetre concernee
if nargin < 2
   hfig = [];
end
if isempty(hfig)
   % cas special
   hfig = get(0,'currentfigure');
   if isempty(hfig)
      return
   end
end

if ~ishandle(hfig)
	return
end
h=getappdata(hfig,'zhandle');
% information pour l'assistant
zuicr(hfig,action);

% detection of cancellation
setappdata(0,'ZGUI_CANCEL',false);

% selon ation
switch lower(action)
	
case {'init','raz'}
%	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	zuireset(h.raz);
	
case {'annulation','close'}
	zuicloseone(hfig);	
	setappdata(0,'ZGUI_CANCEL',true);

case 'validation'
	zuidownloadform(hfig);
	zuisavenonok;
	zuireset(h.validation);
   % code retour
   code_retour = getappdata(hfig,'code_retour');
	% il n'y a pas d'autre variable ici
	zuicloseone(hfig);
   
   if ~isempty(code_retour)
      evalin('base',code_retour);
   end
	
case 'update'
	zuidownloadform(hfig);
	zuisavenonok;
	zuireset(h.update);
   	% code retour
   	code_retour = getappdata(hfig,'code_retour');
   	%disp('update')
   if ~isempty(code_retour)
      evalin('base',code_retour);
   end
   
case 'advanced'
       commande_mem = getappdata(hfig,'commande_mem');
       set(hfig,'visible','off');
       drawnow
       delete(hfig);		
       if isfield(commande_mem,'update')
       		evalin('base',sprintf('zuicreefunform({''%s''},''%s'',%d,%d,''%s'',%d)', ...
		        commande_mem.fonction,commande_mem.racine,commande_mem.nombre,commande_mem.inerte , ...
			commande_mem.code_retour,commande_mem.update));
       else
       		evalin('base',sprintf('zuicreefunform({''%s''},''%s'',%d,%d,''%s'')', ...
		        commande_mem.fonction,commande_mem.racine,commande_mem.nombre,commande_mem.inerte , ...
			commande_mem.code_retour));
       end
	
otherwise
	warning('ation non prise en compte')
	
end

