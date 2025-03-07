% fonction de recherche d'une variable dans l'interface
function  info = zcestou

% dialogue pour selectionner le nom de la varibale
prompt={'field name in the structure (i.e. : ''nequi''), <empty> -> key words :'};
def={''};
dlgTitle='how can we modify this variable ?';
lineNo=1;
answer=zinputdlg(prompt,dlgTitle,lineNo,def);
if isempty(answer)
	return
end

nomvar   = answer{1};
if isempty(nomvar)
	% liste des mots reservee
	s =load('zineb_mots_reserves.mat','reserve');
	reserve = s.reserve;
	clef ={};
	for k = 3:length(reserve)
		if isempty(findstr(reserve{k},'.'))
			if isempty(strmatch(reserve{k},clef,'exact'))
				clef{end+1} = reserve{k};
			end
		end
	end
	clef =sort(clef);
	[s,v] = listdlg('PromptString','Choose a kew word :',...
	              	'SelectionMode','single',...
		            'ListString',clef);
	if v == 0 | isempty(s)
		return
	else
		nomvar =clef{s};
	end
end

% essai de lecture
s =load('zineb_mots_reserves.mat','dico');
dico = s.dico; 
if isfield(dico,nomvar)
	info = zgetfield(dico,nomvar);
	if length(info) >1 
		[s,v] = listdlg('PromptString','Choose a variable name :',...
		'SelectionMode','single',...
		'ListString',info);
		if v == 0 | isempty(s)
			return
		else
			nomvar =info{s};
		end

	elseif iscell(info)
		nomvar =info{1};
	else
		nomvar =info;
	end

	% substitution
	if isempty(findstr(nomvar,'neomulti'))
		nomvar = strrep(nomvar,'param.cons','param.fonction');
	end
	ok = 0;
	root = getappdata(0,'root');
	directory = fullfile(root,'interface','*.m');
	% 1ere recherche
	[s,t] = unix(sprintf('grep %s %s | grep %%',nomvar,directory));
	if s ~= 0 & ~isempty(t)
		warndlg(t,'unix error (grep...)');
		return
	elseif ~ isempty(t)
		% decodage
		tt      = tseparec(t);
		%disp(tt)
		[rep,s] = strtok(t,':');
		if ~isempty(rep) & s(2) == '%'
			[void,rep,ext] = fileparts(rep);
			if ~ isempty(rep)
				ok = 1;
 			end
		end

	end
	if ok == 0
		% deuxieme recherche
		[s,t] = unix(sprintf('grep info.%s %s',nomvar,directory));
		if s ~= 0 & ~isempty(t)
			warndlg(t,'unix error (grep...)');
			return
		elseif ~isempty(t)
			% decodage
			tt      = tseparec(t);
			%disp(tt)
			[rep,s] = strtok(t,':');
			if ~isempty(rep)
				[void,rep,ext] = fileparts(rep);
				if ~isempty(rep)
					ok = 1;
				end
			end
		end
	end
	if ok == 0
		% 3ieme  recherche
		[s,t] = unix(sprintf('grep %s %s',nomvar,directory));
		if s ~= 0 & ~isempty(t)
			warndlg(t,'unix error (grep...)');
			return
		elseif isempty(t)
			warndlg(sprintf('this variable ''%s'' was not found by the search engine.',nomvar),'bad luck !');
			return
		end
		% decodage
		tt      = tseparec(t);
		%disp(tt)
		[rep,s] = strtok(t,':');
		if ~isempty(rep)
			[void,rep,ext] = fileparts(rep);
			if isempty(rep)
				warndlg('No function associated to this variable ...','always unlucky !');
				return
			end
		else
			warndlg('No function associated to this variable ...','you are an unlucky fellow !');
			return
		end
	end
else
	warndlg('undocumented variable ...','game over !');
	return
end



% information
disp('--------------------------------------------------');
info = zinfo;
info = zgetfield(info,nomvar);
fprintf('%s@%s : \n%s \n',rep,nomvar,info);
disp('--------------------------------------------------');

% commande et substitution
if ~isempty(findstr(rep,'zconv'))
	cmd = 'zuiedit_param_gene_funf(''zconv'',''param.gene'',''Convergence'');';
elseif ~isempty(findstr(rep,'zmulti'))
	cmd = 'zuiedit_param_gene_funf(''zmulti'',''param.cons.neomulti'',''Multiplicator'');';
elseif ~isempty(findstr(rep,'zequa'))
	cmd = 'zuiedit_param_gene_funf(''zequa'',''param.gene'',''Equation configuration'');';
elseif ~isempty(findstr(rep,'zexec'))
	cmd = 'zuiedit_param_gene_funf(''zexec'',''param.gene'',''Execution'','''',''zuiedit_param_gene_exec_ctrl'');';
elseif ~isempty(findstr(rep,'zprofile'))
	cmd = 'zuiedit_param_gene_funf(''zprofile'',''param.profile'',''profiler configuration'');';
else
	cmd = sprintf('%s;',rep);
end
if ~isempty(cmd)
	eval(cmd,'1;');
end
