% get list of data for SAL
function tree = get_sal_tree(root_path,authorisation_machine,revision,tree)


if nargin < 1
    root_path = 'data';
    authorisation_machine = 'JET';    
elseif nargin < 2
    authorisation_machine = 'JET';
end
if ~isstruct(authorisation_machine)
    authorisation_machine = getappdata(0,sprintf('SAL_AUTHORISATION_%s',authorisation_machine));
end

% test token validity
if authorisation_machine.validity < now
    error('Token validity has expired: You must renew you authentication');
end

% see details on https://simple-access-layer.github.io/documentation/server.html
if (nargin < 3) || isempty(revision)
    revision = 0; % default provided by the server
end

% initialisation
if nargin < 4
    tree = [];
end


% recursive call off the function on children
data = sal_list(root_path,authorisation_machine,revision);
if isempty(data)
    keyboard
    return
else
    % switch on type
    switch data.type
        case 'branch'
            % loop on children
            for k=1:length(data.object.children)
                if isfield(data.object.children(1),'branches') &&  ...
                   ~isempty(data.object.children.branches)
                        for l=1:length(data.object.children.branches)
                                if all((abs(data.object.children.branches{l}) > 47) & (abs(data.object.children.branches{l}) < 58))
                                    if length(data.object.children) == 1
                                        tree = data;
                                        return
                                    else
                                        tree.(sprintf('children_%d',k)) = data;
                                    end
                                else
                                    loc_name  = strrep(data.object.children.branches{l},'-','_m_');
                                    if (abs(loc_name(1)) > 47) &&  (abs(loc_name(1)) < 58) 
                                        loc_name = sprintf('num_%s',loc_name);
                                    end
                                    tree.(loc_name) =  ...
                                        get_sal_tree(sprintf('%s/%s',root_path,data.object.children.branches{l}), ...
                                        authorisation_machine,revision,tree);
                                end
                        end
                end
                if isfield(data.object.children(1),'leaves') &&  ...
                        ~isempty(data.object.children.leaves)
                    tree = data;
                end
            end
        case 'leaf'
            tree = data;
    end

end






