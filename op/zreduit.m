% ZREDUIT (de)compact ou (de)compress la struture data
%---------------------------------------------------------------------
% fichier zreduit.m ->  zreduit
%
%
% fonction Matlab 5 :
%
% Cette fonction sert a compacter ou a comprimer la strcuture de donnees 
% data. Elle assure aussi les operations inverses. L'opoeration de
% de decompression est lente. 
%  
% syntaxe  :
%  
%     data=zreduit(param,data,mode);
%    
% entree :
%
%     data       =  structure de donnees contenant les variables 
%                   dependant du temps.
%                  
%     param      = structure de donnees contenant les variables 
%                  independantes du temps. 
% 
%     mode       = mode de fonctionnement :
%                     'compact'   -> compact les constantes et remplace les matrices a majorite de 0 par des matrices creuses
%                     'uncompact' -> operation inverse de compact
%                     'compress'  -> compact + suppression des derivee spatiale + suppression des donnees 2D de l'equilibre 
%                     'uncompress' -> operation inverse de compress (il ne faut pas etre presse)
%
% sorties :
% 
%     data       = structure data modifiee
%     
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.0, du 12/02/2003.
% 
% 
% liste des modifications : 
%
% * 09/10/2001 -> suppression des zregular
% * 12/02/2003 -> correction bug sur les complexes
%
%--------------------------------------------------------------
%
function data=zreduit(param,data,mode_reduction)

% mode_reduction = 'compact'   -> supprime les constantes
% mode_reduction = 'uncompact' -> operation inverse de compact
% mode_reduction = 'compress'  -> compact + suppression des derivee spatiale + suppression des donnees 2D de l'equilibre 
% mode_reduction = 'uncompress' -> operation inverse de compress (il ne faut pas etre presse)

% test des entrees
if nargin < 3
	error('Il faut donner une structureparam, une structure data et un mode_reduction')
end
if isempty(data)
	disp('donnees vides : pas d''operations')
	return
end
if isempty(mode_reduction)
	disp('mode_reduction non preciser : pas d''operations')
	return
end


% boucle sur la strcuture
champ = fieldnames(data);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('data.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while (~isempty(champ))
	%premier champ de la liste
	champc=champ{1};
	champ(1)=[];
	eval(strcat('test=isstruct(',champc,');'));
	if test
		% cas d'une sous structure -> ajout de champs
		eval(strcat('champnew=fieldnames(',champc,');'));   	
		for k=1:length(champnew)
			% nom complet pour acceder a la variable
			champnew{k}=strcat(champc,'.',champnew{k});
		end
		% ajout a la liste des champs
		champ=cat(1,champ,champnew);
	else
		% juste pour les tests
		%disp(champc);
		
		% information sur la donnees
		val  = eval(champc,[]); 
		ss   = size(val);
		td   = isa(val,'double');
		tsp  = issparse(val);
		tstr = ischar(val);
		if tstr 
			tstr = ~isempty(findstr(val,'ones('));
		end
		td1  = ~isempty(findstr(champc,'.prof.')) & strcmp(champc((end-1):end),'d1');
		td2  = ~isempty(findstr(champc,'.prof.')) & strcmp(champc((end-1):end),'d2');
		
		% compression/decompression
		if ~isempty(val)
			% selon le mode_reduction
			if (strcmp(mode_reduction,'compact') | strcmp(mode_reduction,'compress')) & (td == 1)
				% si un grand nombre de zeros
				if all(val(:) == val(1)) &  ~iscomplex(val(:))
					val =sprintf('%g .* ones(%s)',val(1),mat2str(ss));
					eval(strcat(champc,' = val;'));
				elseif (length(find(val == 0)) > (0.5 * prod(ss))) & (length(ss) <= 2)
					val = sparse(val);
					eval(strcat(champc,' = val;'));
				end
			elseif (strcmp(mode_reduction,'uncompact') | strcmp(mode_reduction,'uncompress'))
				if tsp == 1		
					val = full(val);
					eval(strcat(champc,' = val;'));
				elseif tstr == 1
					val=eval(val,[]);
					eval(strcat(champc,' = val;'));
				end
			end
			if strcmp(mode_reduction,'compress')
				if (td1 == 1) | (td2 == 1)
					val =sprintf('NaN .* ones(%s)',mat2str(ss));
					eval(strcat(champc,' = val;'));
				end
			end
		end	
	end
end


% autres effets de la compression
if strcmp(mode_reduction,'compress')
	
% 	% seules les donnees de la structure gene non recalculable
% 	gene               = data.gene;
% 	data.gene          = [];
% 	data.gene.temps    = gene.temps;
% 	data.gene.dt       = gene.dt;
% 	data.gene.conv     = gene.conv;
% 	data.gene.nbsplit  = gene.nbsplit;
% 	data.gene.flops    = gene.flops;
% 	data.gene.cputime  = gene.cputime;
% 	data.gene.memory   = gene.memory; 
% 	data.gene.datation = gene.datation;
% 	
	% suppresion des donnees 2D de l'equilibre
	data.equi.R         = sprintf('NaN .* ones(%s)',mat2str(size(data.equi.R)));
	data.equi.Z         = sprintf('NaN .* ones(%s)',mat2str(size(data.equi.Z)));
	data.equi.BR        = sprintf('NaN .* ones(%s)',mat2str(size(data.equi.BR)));
	data.equi.BZ        = sprintf('NaN .* ones(%s)',mat2str(size(data.equi.BZ)));
	data.equi.BPHI      = sprintf('NaN .* ones(%s)',mat2str(size(data.equi.BPHI)));
	data.equi.rhoRZ     = sprintf('NaN .* ones(%s)',mat2str(size(data.equi.rhoRZ)));
	
elseif strcmp(mode_reduction,'uncompress')
	datak =[];
	% boucle sur les temps
	for k=param.gene.kmin:param.gene.k
		% extraction des donnees
		fprintf('Reconstitution de l''indice  %d :\n',k);
		param.gene.t       = data.gene.temps(k);
		param.gene.dt      = data.gene.dt(k);
		if isempty(datak)
			datak              = zget1t(data,k);
		else
			datak              = zget1t(data,k,datak);
		end	
		% appel de l'equilibre
		mem                     = datak.mode.premiertemps;
		ipmem                   = datak.cons.ip;
		datak.mode.premiertemps = 1;
		datak.cons.ip           = datak.equi.ip;
		try
			[cr,equi,void,param.memoire.equi] = ...
                        zequilibre(param.fonction.equi,param.cons.equi,param.memoire.equi,param.gene, ....
                                     param.phys,datak,datak.equi,datak.prof.dpsidt,param.gene.dt,0,1);
		catch
			[cr,equi] = zequilibre(param.fonction.equi,param.cons.equi,param.gene, ....
                                     param.phys,datak,datak.equi,datak.prof.dpsidt,param.gene.dt,0,1);
		end
      datak.mode.premiertemps = mem;
		datak.cons.ip           = ipmem;
		% affectation 
      datak.equi.R         = equi.R;
      datak.equi.Z         = equi.Z;
      datak.equi.BR        = equi.BR;
      datak.equi.BZ        = equi.BZ;
      datak.equi.BPHI      = equi.BPHI;
      datak.equi.rhoRZ     = equi.rhoRZ;
      
      % recalcul des derivees
      disp('recalcul des derivees :')
      datak.prof.psid1     = pdederive(param.gene.x,datak.prof.psi,0,2,2,1);
      datak.prof.psid2     = pdederive(param.gene.x,datak.prof.psi,1,2,2,2);
      
      datak.prof.ned1      = pdederive(param.gene.x,datak.prof.ne,0,2,2,1);
      datak.prof.ned2      = pdederive(param.gene.x,datak.prof.ne,1,2,2,2);
      
      datak.prof.aed1      = pdederive(param.gene.x,datak.prof.ae,0,2,2,1);
      datak.prof.aed2      = pdederive(param.gene.x,datak.prof.ae,1,2,2,2);
      
      datak.prof.ped1      = pdederive(param.gene.x,datak.prof.pe,0,2,2,1);
      datak.prof.ped2      = pdederive(param.gene.x,datak.prof.pe,1,2,2,2);
      
      datak.prof.piond1    = pdederive(param.gene.x,datak.prof.pion,0,2,2,1);
      datak.prof.piond2    = pdederive(param.gene.x,datak.prof.pion,1,2,2,2);
      
      datak.prof.fluced1   = pdederive(param.gene.x,datak.prof.fluce,0,2,2,1);
      datak.prof.fluced2   = pdederive(param.gene.x,datak.prof.fluce,1,2,2,2);
      
      datak.prof.fluciond1 = pdederive(param.gene.x,datak.prof.flucion,0,2,2,1);
      datak.prof.fluciond2 = pdederive(param.gene.x,datak.prof.flucion,1,2,2,2);
      
      datak.prof.rotd1   = pdederive(param.gene.x,datak.prof.rot,0,2,2,1);
      datak.prof.rotd2   = pdederive(param.gene.x,datak.prof.rot,1,2,2,2);
      
%       datak.prof.ned1      = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.ne,0,2,2,1));
%       datak.prof.ned2      = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.ne,1,2,2,2));
%       
%       datak.prof.aed1      = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.ae,0,2,2,1));
%       datak.prof.aed2      = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.ae,1,2,2,2));
%       
%       datak.prof.ped1      = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.pe,0,2,2,1));
%       datak.prof.ped2      = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.pe,1,2,2,2));
%       
%       datak.prof.piond1    = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.pion,0,2,2,1));
%       datak.prof.piond2    = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.pion,1,2,2,2));
%       
%       datak.prof.fluced1   = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.fluce,0,2,2,1));
%       datak.prof.fluced2   = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.fluce,1,2,2,2));
%       
%       datak.prof.fluciond1 = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.flucion,0,2,2,1));
%       datak.prof.fluciond2 = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.flucion,1,2,2,2));
%       
%       datak.prof.rotd1   = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.rot,0,2,2,1));
%       datak.prof.rotd2   = zregular(param.gene.x,pdederive(param.gene.x,datak.prof.rot,1,2,2,2));
      
      % on recolle les morceaux
      data = zput1t(data,k,datak);
      
	end	
end
