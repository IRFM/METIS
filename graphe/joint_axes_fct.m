function joint_axes_fct(action)

% gestion des entrees
if nargin < 1
    return
elseif isempty(action)
    return
end

% constantes
h = gcf;


switch action
    
    case 'xmin'
        [x,y] = ginput(1);
        hh    = findobj(h,'type','axes');
        for k =1:length(hh)
            xlim = get(hh(k),'xlim');
            if x < xlim(2)
                xlim(1) = x;
            end
            hl = legend(hh(k));
            if isempty(hl) 
                    set(hh(k),'xlim',xlim);
            elseif strcmp(get(hh(k),'tag'),'legend')
                %rien
            else
                legend(hh(k),'hide');
                set(hh(k),'xlim',xlim);
                legend(hh(k),'show');
            end
        end
    case 'xmax'
        [x,y] = ginput(1);
        hh    = findobj(h,'type','axes');
        for k =1:length(hh)
            xlim = get(hh(k),'xlim');
            if x > xlim(1)
                xlim(2) = x;
            end
            hl = legend(hh(k));
            if isempty(hl) 
                    set(hh(k),'xlim',xlim);
            elseif strcmp(get(hh(k),'tag'),'legend')
                %rien
            else
                legend(hh(k),'hide');
                set(hh(k),'xlim',xlim);
                legend(hh(k),'show');
            end
        end
     case 'auto'
        hh    = findobj(h,'type','axes');
        for k =1:length(hh)
            hl = legend(hh(k));
            if isempty(hl) 
                  set(hh(k),'xlimmode','auto','ylimmode','auto');
            elseif strcmp(get(hh(k),'tag'),'legend')
                %rien
            else
                legend(hh(k),'hide');
                set(hh(k),'xlimmode','auto','ylimmode','auto');
                legend(hh(k),'show');
            end
            hl = legend(hh(k));
            if isempty(hl)
                 set(hh(k),'xlimmode','auto','ylimmode','auto');
            end
        end        
    otherwise
        disp('Action non prise en compte')
end
zoom yon


