% ZUIPARAM_FCT callback du formulaire de parametrage des profils
%---------------------------------------------------------------
% fichier zuiparam_fct.m 
%
% fonction Matlab 5 :
% Cette fonction gere les "callback" du formulaire de parametrage des profils
%              zuiparam.m
%  
% syntaxe  :
%	zuiparam_fct(action);
%    
% entree :
%	action = chaine de caratere, nom de l'action
%
% sorties : aucune
% 
% fonction ecrite par C. Passeron, psote 61 19
% version 1.6, du 11/09/2001.
% 
% liste des modifications : 
%   * XX/XX/20XX -> commentaires  
%
%--------------------------------------------------------------
%
function zuiparam_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hform,hui] = zuiformhandle('parametrer') ;

x      = evalin('base','param.gene.x') ;
centre = zuidata(hui.edit_centre) ;
bord   = zuidata(hui.edit_bord) ;
alpha  = zuidata(hui.edit_alpha) ;
beta   = zuidata(hui.edit_beta) ;

% information pour l'assistant
zuicr(hform,action) ;

% selon ation
switch lower(action)
	
case 'edit_centre'
	
case 'edit_bord'

case 'edit_alpha'
	if alpha<=0
		st = sprintf(' ALPHA must be > 0') ;
		disp(st)
		defaut = getappdata(hui.edit_alpha,'init_defaut') ;
		zuidata(hui.edit_alpha,defaut.string) ;
		return
	end

case 'edit_beta'
	if beta<=0
		st = sprintf(' BETA must be > 0') ;
		disp(st)
		defaut = getappdata(hui.edit_beta,'init_defaut') ;
		zuidata(hui.edit_beta,defaut.string) ;
		return
	end
		
case {'init','raz'}
	zuiformvisible(hform) ;
	zuiformreset(hform) ;
	zuiuploadform(hform) ;
	centre = zuidata(hui.edit_centre) ;
	bord   = zuidata(hui.edit_bord) ;
	alpha  = zuidata(hui.edit_alpha) ;
	beta   = zuidata(hui.edit_beta) ;
	zuireset(hui.raz)
	
case {'annulation','close'}
	variable = getappdata(hform,'variable') ;
	switch variable
	case 'data.prof.te'
		evalin('base','zrecalcultemppres(''el'',''clear'');');
	case 'data.prof.pe'
		evalin('base','zrecalcultemppres(''el'',''clear'');');
	case 'data.prof.ti'
		evalin('base','zrecalcultemppres(''ion'',''clear'');');
	case 'data.prof.pion'
		evalin('base','zrecalcultemppres(''ion'',''clear'');');
	end	
	zuireset(hui.annulation)
	zuicloseone(hform) ;
	return
	
case 'validation'
	% recupration de la modulation
	modulation = get(hui.pop_modul,'value') ;
	var_modul  = getappdata(hform,'var_modul') ;
	tm     = evalin('base',var_modul) ;
	
	% suppression des NaN et inf
	ind    = find(~isfinite(tm));
	if ~isempty(ind)
		tm(ind) = zeros(1,length(ind));
	end
	% sommation du  module sur l'espace/voies
	tm     = sum(abs(tm),2);
			
	% normalisation du max a 1
	tmax   = max(tm);
	if tmax == 0
		tm = ones(size(tm));
	else
		tm = tm ./ tmax;
	end
	
	profil = 	(centre - bord) * (1 - x.^alpha).^beta + bord ;

	% selon le cas constant ou module
	if modulation == 1
	           	profil = tm * profil ; 
	else
	           	profil = ones(size(tm)) * profil ;
	end

	variable = getappdata(hform,'variable') ;
	nbt = evalin('base','param.gene.nbt') ;
	zassignin('base',variable,profil) ;
	switch variable
	case 'data.prof.te'
		evalin('base','zrecalcultemppres(''el'');');
	case 'data.prof.pe'
		evalin('base','zrecalcultemppres(''el'');');
	case 'data.prof.ti'
		evalin('base','zrecalcultemppres(''ion'');');
	case 'data.prof.pion'
		evalin('base','zrecalcultemppres(''ion'');');
	end	
	zuireset(hui.validation) ;
	zuicloseone(hform) ;
	zuisavenonok;
	return

otherwise
%      warning('ation non prise en compte')
end

% calcul du nouveau profil 
profil = (centre - bord) * (1 - x.^alpha).^beta + bord ;

% trace du profil
axes(hui.axes_plot)
plot(x,profil)
title('(center-edge)*(1-x^\alpha)^\beta + edge')
xlabel('x (su)');
zoom on

