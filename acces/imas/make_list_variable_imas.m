% cette fonction deplace le temps de la racine de la structure matlab vers les champs de donnees
% s_in est une structure compatible avec l'UAL (s_in(temps).data(espace))
% s_out est une structure de type cronos (s_out.data(temps,espace)) 
function make_list_variable_imas(nom_root)

% generate IDS
s_in = ids_gen(nom_root);

% recherche du nombre de temps
nbt = length(s_in);

% boucle sur les temps
for lk = 1:nbt
    if iscell(s_in)
        data = s_in{lk};
    else
        data = s_in(lk);
    end
    
    % liste des champs la structures
    champ = fieldnames(data);
    for k=1:length(champ)
        % nom complet pour acceder a la variable
        champ{k}=sprintf('data.%s',champ{k});
    end
    
    % jusqu'a ce qu'il n'y ait plus de champ
    test=0;
    while(~isempty(champ))
        %premier champ de la liste
        champc=champ{1};
        champ(1)=[];
        locvar = [];
        locvar = eval(champc);
        %eval(sprintf('test=isstruct(%s);',champc));
        test = isstruct(locvar);
        test_cell = iscell(locvar);
        len  = length(locvar);        
        if test && (len == 1)
            % cas d'une sous structure -> ajout de champs
            champnew=fieldnames(locvar);
            %eval(sprintf('champnew=fieldnames(%s);',champc));
            for k=1:length(champnew)
                % nom complet pour acceder a la variable
                champnew{k}=sprintf('%s.%s',champc,champnew{k});
            end
            % ajout a la liste des champs
            if isempty(champ)
                champ =champnew;
            else
                champ=cat(1,champ,champnew);
            end
        elseif  test_cell && (len == 1)
            % cas d'une sous structure -> ajout de champs
            champnew=fieldnames(locvar{1});
            %eval(sprintf('champnew=fieldnames(%s);',champc));
            for k=1:length(champnew)
                % nom complet pour acceder a la variable
                champnew{k}=sprintf('%s{1}.%s',champc,champnew{k});
            end
            % ajout a la liste des champs
            if isempty(champ)
                champ =champnew;
            else
                champ=cat(1,champ,champnew);
            end
        elseif test && (len > 0)
            % cas autre de structure indexee
            champloc = fieldnames(locvar);
            champnew = {};
            for k=1:length(champloc)
                for kzz  = 1:len
                    % nom complet pour acceder a la variable
                    champnew{end+1}=sprintf('%s(%d).%s',champc,kzz,champloc{k});
                end
            end
            % ajout a la liste des champs
            if isempty(champ)
                champ =champnew';
            else
                champ=cat(1,champ,champnew');
            end
        elseif test_cell && (len > 0)
            % cas autre de structure indexee
            for lz =1:length(locvar)
                if isstruct(locvar{lz})
                    champloc = fieldnames(locvar{lz});
                    break;
                end
            end
            champnew = {};
            for k=1:length(champloc)
                for kzz  = 1:len
                    % nom complet pour acceder a la variable
                    champnew{end+1}=sprintf('%s{%d}.%s',champc,kzz,champloc{k});
                end
            end
            % ajout a la liste des champs
            if isempty(champ)
                champ =champnew';
            else
                champ=cat(1,champ,champnew');
            end
           
            
        elseif isempty(findstr(champc,'_error_'))
            
            % juste pour les tests
           disp(strrep(champc,'data.',sprintf('%s.',nom_root)));
            
        end
    end
end
