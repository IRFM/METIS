%  Z0LOGLIN  courte description  
%------------------------------------------------------------------------------- 
% fichier :  z0loglin.m  ->  z0loglin 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   z0loglin(ha) 
%  
% entrees :  
%  ha = 
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
function z0loglin(ha)

if nargin >0
	pos  = get(ha,'pos');
	x    = pos(1);
	y    = pos(2) + pos(4);
	dx   = max(0.001,pos(4)/7);
	dy   = dx ./ 2;
	posc = [x,y-dy,dx,dy];
	uicontrol(get(ha,'parent'),'style','togglebutton','callback','z0loglin', ...
	          'units',get(ha,'units'),'position',posc,'userdata',ha);

else
	ho = gco;
	ha = get(ho,'userdata');
	sc = get(ha,'yscale');
	switch sc
	case 'linear'
		set(ha,'yscale','log');
	otherwise
		set(ha,'yscale','linear');
	end
end
