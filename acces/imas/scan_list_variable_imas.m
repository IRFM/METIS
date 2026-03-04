function scan_list_variable_imas(s_in,info_in,nom_root,fid)

% recherche du nombre de temps
nbt = length(s_in);
if nbt > 10
    nbt = 1;
end
% boucle sur les temps
for lk = 1:nbt
    if iscell(s_in)
        data = s_in{lk};
        info = info_in;
    else
        data = s_in(lk);
        info = info_in;
    end
    
    % liste des champs la structures
    champ      = fieldnames(data);
    champ_info = fieldnames(info);
    for k=1:length(champ)
        % nom complet pour acceder a la variable
        champ_info{k}=sprintf('info.%s',champ{k});
        champ{k}=sprintf('data.%s',champ{k});
    end
    
    % jusqu'a ce qu'il n'y ait plus de champ
    test=0;
    while(~isempty(champ))
        %premier champ de la liste
        champc=champ{1};
        champ(1)=[];
        if isempty(champ_info) || ~iscell(champ_info)
            continue;
        end
        champ_infoc = champ_info{1};
        champ_info(1)=[];
        locvar = [];
        locinfo = '';
        try
            locvar = eval(champc);
            try
                locinfo = eval(champ_infoc);
                %eval(sprintf('test=isstruct(%s);',champc));
                test = isstruct(locvar);
                test_cell = iscell(locvar);
                len  = length(locvar);
                if test && (len == 1)
                    % cas d'une sous structure -> ajout de champs
                    if ~isempty(locinfo)
                        champnew=fieldnames(locvar);
                        try
                            champ_infonew=fieldnames(locinfo);
                            %eval(sprintf('champnew=fieldnames(%s);',champc));
                            for k=1:length(champnew)
                                % nom complet pour acceder a la variable
                                champ_infonew{k}=sprintf('%s.%s',champ_infoc,champnew{k});
                                champnew{k}=sprintf('%s.%s',champc,champnew{k});
                            end
                            % ajout a la liste des champs
                            if isempty(champ)
                                champ =champnew;
                                champ_info = champ_infonew;
                            else
                                champ=cat(1,champ,champnew);
                                champ_info=cat(1,champ_info,champ_infonew);
                            end
                        catch
                            % this no more a data structure in current DD
                        end
                    end
                    
                elseif  test_cell && (len == 1)
                    % cas d'une sous structure -> ajout de champs
                    if ~isempty(locinfo)
                        champnew=fieldnames(locvar{1});
                        champ_infonew=fieldnames(locinfo);
                        %eval(sprintf('champnew=fieldnames(%s);',champc));
                        for k=1:length(champnew)
                            % nom complet pour acceder a la variable
                            champ_infonew{k}=sprintf('%s.%s',champ_infoc,champnew{k});
                            champnew{k}=sprintf('%s{1}.%s',champc,champnew{k});
                        end
                        % ajout a la liste des champs
                        if isempty(champ)
                            champ =champnew;
                            champ_info = champ_infonew;
                        else
                            champ=cat(1,champ,champnew);
                            champ_info=cat(1,champ_info,champ_infonew);
                        end
                    end
                elseif test && (len > 0)
                    % cas autre de structure indexee
                    champloc = fieldnames(locvar);
                    champnew = {};
                    champ_infonew={};
                    for k=1:length(champloc)
                        for kzz  = 1:len
                            % nom complet pour acceder a la variable
                            champnew{end+1}=sprintf('%s(%d).%s',champc,kzz,champloc{k});
                            champ_infonew{end+1}=sprintf('%s.%s',champ_infoc,champloc{k});
                        end
                    end
                    % ajout a la liste des champs
                    if isempty(champ)
                        champ =champnew';
                        champ_info = champ_infonew';
                    else
                        champ=cat(1,champ,champnew');
                        champ_info=cat(1,champ_info,champ_infonew');
                    end
                elseif test_cell && (len > 0)
                    % not all the input if too long
                    len = 1;
                    % cas autre de structure indexee
                    for lz =1:length(locvar)
                        if isstruct(locvar{lz})
                            champloc = fieldnames(locvar{lz});
                            break;
                        end
                    end
                    champnew = {};
                    champ_infonew={};
                    for k=1:length(champloc)
                        for kzz  = 1:len
                            % nom complet pour acceder a la variable
                            champnew{end+1}=sprintf('%s{%d}.%s',champc,kzz,champloc{k});
                            champ_infonew{end+1}=sprintf('%s.%s',champ_infoc,champloc{k});
                        end
                    end
                    % ajout a la liste des champs
                    if isempty(champ)
                        champ =champnew';
                        champ_info = champ_infonew';
                    else
                        champ=cat(1,champ,champnew');
                        champ_info=cat(1,champ_info,champ_infonew');
                    end
                    
                    
                else
                    % test if field is not empty or NaN
                    try
                        rep = eval(champc);
                    catch
                        % Not a IMAS data in present version of DD
                        rep = [];
                    end
                    isfilled = false;
                    if ~isempty(rep) && any(isfinite(rep(:)))
                        if isnumeric(rep)
                            if all(rep == fix(rep))
                                if any(rep(:) > -999999999)
                                    isfilled = true;
                                end
                            elseif any(rep(:) > -9.0000e+40)
                                isfilled = true;
                            end
                        else
                            isfilled = true;
                        end
                    end
                    if isfilled
                        try
                            infodisp = eval(champ_infoc);
                        catch
                            infodisp = '';
                        end
                        % test if it is a METIS IDS data from previous IMAS version
                        if ~isempty(infodisp)
                            % juste pour les tests
                            fprintf(fid,'* %s: %s\n',strrep(strrep(strrep(champc,'data.',sprintf('%s.',nom_root)),'{1}','{:}'),'.','/'),infodisp);
                        end
                    end
                    
                end
            catch
                locinfo = '';
            end
            
        catch
            % not a present version IDS data
        end
    end
    
end