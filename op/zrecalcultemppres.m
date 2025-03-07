%  ZRECALCULTEMPPRES  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zrecalcultemppres.m  ->  zrecalcultemppres 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   zrecalcultemppres(mode,clear) 
%  
% entrees :  
%  mode  = 
%  clear = 
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
function zrecalcultemppres(mode,clear)

physe = evalin('base','param.phys.e','[]');
if isempty(physe)
	disp('Qu''avez vous fait des donnees ?')
	return
end

switch mode
case 'el'
 	if nargin > 1
		evalin('base','param.edit.tepe ='''';');
		return
	end
	tepe = evalin('base','param.edit.tepe');
	dens = evalin('base','data.prof.ne');
	if all(~isfinite(dens))
		disp('La densite electronique n''est pas definie => pas de recalcul de Te ou Pe');
		return
	end
   switch tepe
	case 'te'
		% reclacul de pe
		pe = evalin('base','data.prof.te') .* dens .* physe;
		zassignin('base','data.prof.pe',pe);
		evalin('base','param.edit.tepe ='''';');
		
	case 'pe'
		% reclacul de te
		te = evalin('base','data.prof.pe') ./ dens ./ physe;
		zassignin('base','data.prof.te',te);
		evalin('base','param.edit.tepe ='''';');
	otherwise
		error('C''est nouveau ? -> tepe non defini');
	end	

case 'ion'
 	if nargin > 1
		evalin('base','param.edit.tipion ='''';');
		return
	end
	tipion = evalin('base','param.edit.tipion');
	dens = evalin('base','data.prof.ni');
	if all(~isfinite(dens))
		dens = evalin('base','data.prof.ne') .*  evalin('base','data.prof.ae');
		if all(~isfinite(dens))
			disp('La densite ionique n''est pas definie => pas de recalcul de Ti ou Pion');
			return
		end
	end
   switch tipion
	case 'ti'
		% reclacul de pion
		pion = evalin('base','data.prof.ti') .* dens .* physe;
		zassignin('base','data.prof.pion',pion);
		evalin('base','param.edit.tipion ='''';');
		
	case 'pion'
		% reclacul de ti
		ti = evalin('base','data.prof.pion') ./ dens ./ physe;
		zassignin('base','data.prof.ti',ti);
		evalin('base','param.edit.tipion='''';');
	otherwise
		error('C''est nouveau ? -> tipion non defini');
	end	

otherwise
	error('Ca ne doit pas arriver -> mode non defini');
end
