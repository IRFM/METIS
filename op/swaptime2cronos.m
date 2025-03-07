% cette fonction deplace le temps de la racine de la structure matlab vers les champs de donnees
% s_in est une structure compatible avec l'UAL (s_in(temps).data(espace))
% s_out est une structure de type cronos (s_out.data(temps,espace)) 
function s_out = swaptime2cronos(s_in)

% recherche du nombre de temps
nbt = length(s_in);
s_out = [];

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
        try
            locvar = eval(champc);
        catch
            % for ggd data
            locvar = [];
            last_p = max(findstr(champc,'.'));
            racine = champc(1:(last_p - 1));
            reste  = champc((last_p + 1):end);
            len_racine = eval(sprintf('length(%s)',racine));
            if eval(sprintf('~isfield(%s,''%s'')',racine,reste))
                % rien
            elseif eval(sprintf('iscell(%s)',racine))
                for ko = 1:len_racine
                    locvar(ko) = eval(sprintf('%s{%d}.%s',racine,ko,reste));
                end    
            else
                for ko = 1:len_racine
                    locvar(ko) = eval(sprintf('%s(%d).%s',racine,ko,reste));
                end                   
            end
        end    
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
           
            
        else
            
            % juste pour les tests
            % disp(champc);
            
            % chaines utiles
            sn = strrep(champc,'data.','s_out.');
            %if lk == 1
            %	eval(sprintf('%s = [];',sn));
            %end
            % creation du champs dans la structure de sorttie
            %eval(sprintf('vide = isempty(%s);',champc));
            %eval(sprintf('isanum = isnumeric(%s);',champc));
            %eval(sprintf('ss_in = size(%s);',champc));
            if ~isempty(locvar)
                if isnumeric(locvar);
                    %locvar = squeeze(locvar);
                    ss_in = size(locvar);
                    nbs = sum((ss_in > 1));
                    if lk == 1
                        if nbs < 2
                            eval(sprintf('%s  = NaN .* ones(cat(2,nbt,length(locvar)));',sn));
                        else
                            eval(sprintf('%s  = NaN .* ones(cat(2,nbt,size(locvar)));',sn));
                        end
                    end
                    switch nbs
                        case 0
                            eval(sprintf('%s(lk,:) = locvar;',sn));
                        case 1
                            eval(sprintf('sold_ = size(%s);',sn));
                            if sold_(2) ~= length(locvar)
                                eval(sprintf('toberedim_ = %s;',sn));
                                newsize_ = max(sold_(2),length(locvar));
                                newvar_  = NaN * ones(sold_(1),newsize_);
                                newvar_(:,1:sold_(2)) = toberedim_;
                                eval(sprintf('%s = newvar_;',sn));
                                newvar_ = NaN * ones(1,newsize_);
                                newvar_(1:length(locvar)) = locvar;
                                locvar = newvar_;
                            end
                            if size(locvar,1) > 1
                                eval(sprintf('%s(lk,:) = locvar;',sn));
                            else
                                eval(sprintf('%s(lk,:) = locvar'';',sn));
                            end
                        case 2
                            repval = shiftdim(locvar,-1);
                            eval(sprintf('%s(lk,:,:) = repval;',sn));
                        case 3
                            repval = shiftdim(locvar,-1);
                            eval(sprintf('%s(lk,:,:,:) = repval;',sn));
                        case 4
                            repval = shiftdim(locvar,-1);
                            eval(sprintf('%s(lk,:,:,:,:) = repval;',sn));
                        case 5
                            repval = shiftdim(locvar,-1);
                            eval(sprintf('%s(lk,:,:,:,:,:) = repval;',sn));
                        case 6
                            repval = shiftdim(locvar,-1);
                            eval(sprintf('%s(lk,:,:,:,:,:,:) = repval;',sn));
                        otherwise
                            error('this matrix dimension is not yet implemented');
                    end
                elseif ischar(locvar)
                    eval(sprintf('%s(lk) = {locvar};',sn));
                    
                elseif lk == 1
                    %eval(sprintf('%s = %s;',sn,champc));
                    eval(sprintf('%s = locvar;',sn));
                end
            end
        end
    end
end
