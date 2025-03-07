%ZINEBALLTEST  test de tous les module de CRONOS 
%------------------------------------------------------------------------------- 
% fichier :  zineballtest.m  ->   zineballtest 
% 
% 
% fonction Matlab 7 : 
% 
% Cette fontction cree un test ou execute un test pour tous les modules de CRONOS 
%  
% syntaxe :  
%
%     creation du test :
%
%   		test = zineballtest(data,param,temps,{type});
%
%     execution du test:
%
%   		test = zineballtest(test);
%  
% entrees :  
%  		data  = structure data de CRONOS
%               param = structure param de CRONOS
%               temps = temps choisi pour le test (s)
%               type  = argument optionnel donnant le type de test
%
%               test  = structure test creee par zineballtest
%    
%    type de test :
%               all     : cree des test pour tous les modules et le solver de CRONOS
%               modules : cree seulement les tests des modules externes
% 
%  
% sorties :  
%  
%		test = structure test creee par zineballtest 
%
% fonction ecrite par J-F Artaud , poste 62-15  
% version  3.1  du  25/10/06  
%  
%  
% liste des modifications :  version CVS
%-------------------------------------------------------------------------------  
%  
function  [test,tolerance] = zineballtest(varargin)

% initialisation
test = [];

% selon le nombre d'entree
if length(varargin) == 1
	mode = 'test';
else
	mode = 'make';
end

switch mode
case 'make'
	% sortie
	test = [];
	% cree la structure test
	data  = varargin{1};
	param = varargin{2};
	temps = varargin{3};
	if length(varargin) >= 4
		type  = varargin{4};
	else
		% valeur possible de type : modules,all
		type  = 'all';
	end
	
	% lecture de la liste des modules	
	[liste_all,liste_module] = zineb1test;
	% choix de la liste
	switch type
	case 'all'
		liste = liste_all;
		
	case 'modules'
		liste = liste_module;
	
	otherwise
		error('unknown list of modules/actions');
	end
	
	% boucle sur la liste
	namelist = fieldnames(liste);
	for k=1:length(namelist)
		module = namelist{k};
		[data_out,param_out,tolerance] = zineb1test(data,param,temps,module);
		if ~isempty(data_out)
			test = setfield(test,module,'data',data_out);
			test = setfield(test,module,'param',param_out);
			test = setfield(test,module,'tolerance',tolerance);
		end
	end

case 'test'

	% boucle sur les test
	test  = varargin{1};
	namelist = fieldnames(test);
	for k=1:length(namelist)
		module = namelist{k};
		data   = getfield(test,module,'data');
		param  = getfield(test,module,'param');
		if isfield(getfield(test,module),'tolerance')
			tolerance = getfield(test,module,'tolerance');
		else
			tolerance = 1e-3;
		end
		[data_out,param_out,tolerance] = zineb1test(data,param,[],module,'tolerance',tolerance);
		test = setfield(test,module,'data',data_out);
		test = setfield(test,module,'param',param_out);
		test = setfield(test,module,'tolerance',tolerance);
	end


otherwise
	disp('????????')
end


