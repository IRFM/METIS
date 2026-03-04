function out = decode_codeparam(codeparam)

% initialisation
out = [];
varm1 = '';
varm2 = '';
% separation
rep = split(codeparam,{newline,'<','>'});
% loop on cell
for k=1:length(rep)
    var = strtrim(rep{k});
    if isempty(var)
        %skip
    else
        if var(1) == '/'
           nom = var(2:end);
           if strcmp(nom,varm2)
               if strcmpi(varm1,'NAN')
                   out.(nom) = NaN;
               elseif strcmpi(varm1,'INF')
                   out.(nom) = Inf;
               else
                   out.(nom) = str2double(varm1);
                   if isempty(out.(nom)) || ~isfinite(out.(nom))
                        out.(nom) = varm1;
                   end
               end
           end
        elseif var(end) == '/'
            % empty variable
            out.(var(1:end-1)) = '';
        end
    end
    varm2 = varm1;
    varm1 = var;
end

