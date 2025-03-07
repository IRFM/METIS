% cette fonction retire les donnees externe de METIS
function  rm_cs4m(menu_mode)

if nargin == 0
    noms = fieldnames(getappdata(0));
    for k = 1:length(noms)
        if findstr(noms{k},'_EXP')
            rmappdata(0,noms{k});
            fprintf('removing external data in METIS for %s\n',strtok(noms{k},'_'));
        end
    end
else
    % create menu    
    noms = fieldnames(getappdata(0));
    menu_item = {};
    for k = 1:length(noms)
        if findstr(noms{k},'_EXP')
            indcut = max(find(noms{k} == '_'));
            menu_item{end+1} = noms{k}(1:indcut-1);
        end
    end
    if isempty(menu_item)
        return
    end
    menu_item{end+1} = 'All';
    menu_item{end+1} = 'Cancel';
    k = menu('Which external data do-you want remove ?',menu_item);
    if k == 0
         return
    elseif k == length(menu_item)
         return
    elseif k  == (length(menu_item) - 1)
        rm_cs4m;
    else
         rmappdata(0,sprintf('%s_EXP',menu_item{k}));
         fprintf('removing external data in METIS for %s\n',menu_item{k});
    end
    
end



