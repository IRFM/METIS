% ZUIEDIT_MODE_PLOT   trace du mode
%--------------------------------------------------------------
% fichier zuiedit_mode_plot.m  
%
% fonction Matlab 5 :
%	trace du mode
%
% syntaxe :
%	zuiedit_mode_plot(x,y)
% 
% entrees :
%  x = valeurs des abscisses  
%  y = valeurs du mode
%
% sorties
%
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 12/07/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
function zuiedit_mode_plot(x,y)

% recupere le handle de la fenetre concernee
%[hfig,h] = zuiformhandle('mode') ;
hfig = gcf;
h    = getappdata(hfig,'zhandle');

% on recupere les valeurs possibles du mode.
valeur        = getappdata(hfig,'valeur') ;

ydata = y ;

for i=1:length(valeur)
	ind = find(y~=valeur{i}) ;
	ydata(ind) = NaN ;

	ind = find(y==valeur{i}) ;
	ydata(ind) = y(ind) ;

	set(getfield(h,strcat('line_',num2str(i))),'xdata',x,'ydata',ydata, ...
		     'visible','on','Marker','o','MarkerSize',3) ;	

end


