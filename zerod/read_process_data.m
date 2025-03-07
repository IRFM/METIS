function out = read_process_data(filename)

out.filename = filename;
fid = fopen(filename);
if fid < 0
    fprintf('uanble to open file %s\n',filename);
    return
end

% boucle de lecture
tk = 1;
no = 1;
section = 'header';
for k=1:100000
    line = fgetl(fid);
    if ~ischar(line)
        break
    elseif ~isempty(strfind(line,'End of PROCESS Output'))
        break
    end
    if (length(strfind(line,'"')) == 2) && isempty(strfind(line,','))
        section = strtrim(strrep(line,'"',''));
    elseif ~isempty(strfind(line,'* Numerics *'))
        section = 'numerics';
    elseif ~isempty(strfind(line,'* Power Reactor Costs *'))
        section = 'power_reactor_cost';
    elseif ~isempty(strfind(line,'* Plant Availability *'))
        section = 'plant_availability';
    elseif ~isempty(strfind(line,'* Plasma *'))
        section = 'plasma';
    elseif ~isempty(strfind(line,'* Current Drive System *')) || ~isempty(strfind(line,'* Current Drive System *'))
        section = 'hcd';
    elseif ~isempty(strfind(line,'* Pulsed Reactor *'))
        section = 'pulsed';
    elseif ~isempty(strfind(line,'* Times *'))
        section = 'times';
    elseif ~isempty(strfind(line,'* Divertor *'))
        section = 'divertor';
    elseif ~isempty(strfind(line,'* Radial Build *'))
        section = 'radial_build';
    elseif ~isempty(strfind(line,'* Vertical Build *'))
        section = 'vertical_build';
    elseif ~isempty(strfind(line,'* TF Coils *'))
        section = 'tf_coils';
    elseif ~isempty(strfind(line,'* Superconducting TF Coils *'))
        section = 'sc_tf_coils';
    elseif ~isempty(strfind(line,'* Central Solenoid and PF Coils *'))
        section = 'cs_pf';
    elseif ~isempty(strfind(line,'* Volt Second Consumption *'))
        section = 'volt_second';
    elseif ~isempty(strfind(line,'* Waveforms *'))
        section = 'waveforms';
    elseif ~isempty(strfind(line,'* Support Structure *'))
        section = 'support_structure';
    elseif ~isempty(strfind(line,'* PF Coil Inductances *'))
        section = 'inductances';
    elseif ~isempty(strfind(line,'* Shield / Blanket *'))
        section = 'shield_blanket';
    elseif ~isempty(strfind(line,'* Superconducting TF Coil Power Conversion *'))
        section = 'tf_power_supplies';
    elseif ~isempty(strfind(line,'* PF Coil Power Conversion *'))
        section = 'pf_power_supplies';
    elseif ~isempty(strfind(line,'* Vacuum System *'))
        section = 'vacuum_system';
    elseif ~isempty(strfind(line,'* Plant Buildings System *'))
        section = 'plant_buildings_system';
    elseif ~isempty(strfind(line,'* AC Power *'))
        section = 'ac_power';
    elseif ~isempty(strfind(line,'* Plant Power / Heat Transport Balance *'))
        section = 'bop';
    elseif  (~isempty(strfind(line,' (')) || ~isempty(strfind(line,'fimp('))) && ~isempty(strfind(line,') '))
        [var_name,value,comment,unit,no] = decode_line_var(line,no);
        try
            out.(section).(var_name).value   = value;
            out.(section).(var_name).comment = comment;
            out.(section).(var_name).unit    = unit;
            out.(section).(var_name).line    = line;
        catch
            [var_name,value,comment,unit,no] = decode_line_var_alt(line,no);
            out.(section).(var_name).value   = value;
            out.(section).(var_name).comment = comment;
            out.(section).(var_name).unit    = unit;
            out.(section).(var_name).line    = line;
        end
    else
        [var_name,value,comment,unit,no] = decode_line_var_alt(line,no);
        if ~isempty(var_name) && ~isempty(value) && (var_name(1) == '_')
            var_name = sprintf('US_%s',var_name(2:end));
        end
        if ~isempty(var_name) 
            var_name = strrep(var_name,'{','_');
            var_name = strrep(var_name,'\','_');
            var_name = strrep(var_name,'}','_');
        end        
        if ~isempty(var_name) && ~isempty(value)
            if ~isfield(out,section)
                out.(section).(var_name).value   = value;
            elseif isfield(out.(section),var_name)
                if iscell(out.(section).(var_name).value)
                    out.(section).(var_name).value{end+1} = value;
                elseif size(out.(section).(var_name).value,1)   == size(value,1)
                    out.(section).(var_name).value = cat(2,out.(section).(var_name).value,value);
                elseif size(out.(section).(var_name).value,2)   == size(value,2)
                    out.(section).(var_name).value = cat(1,out.(section).(var_name).value,value);
                else
                    out.(section).(var_name).value{1} = out.(section).(var_name).value{1};
                    out.(section).(var_name).value{end+1} = value;
                end
            else
                out.(section).(var_name).value   = value;
            end
            out.(section).(var_name).comment = comment;
            out.(section).(var_name).unit    = unit;
            out.(section).(var_name).line    = line;
        elseif isfield(out,section)
            out.(section).(sprintf('line_%d',length(fieldnames(out.(section)))+1)) = line;
        else
            out.(section).(sprintf('line_%d',1)) = line;
        end
    end
end

fclose(fid);


function [name,value,comment,unit,no] = decode_line_var(line,no)


ind_deb = strfind(line,' (');
if isempty(ind_deb)
    ind_deb = strfind(line,'fimp(')-2;
    unit_type = 1;
else
    unit_type = 0;
end
ind_fin = strfind(line,') ');

if unit_type == 1
    unit = line(1:3);
    unit(unit ==' ') = '';
    unit(unit =='_') = '';
elseif (length(ind_deb) > 1)  && (length(ind_fin) > 1)
    unit = line((ind_deb(1)+2):(ind_fin(1) -1));
else
    unit = '';
end
name = line((ind_deb(end)+2):(ind_fin(end) -1));
name(name <= ' ') = '';
indr= find(name == '=');
if ~isempty(indr)
    name(indr) = '';
end
indr= find(name == '%');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == ',');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == '-');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == '+');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == '*');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == '.');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == ':');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == '/');
if ~isempty(indr)
    name(indr) = 'o';
end
indr= find(name == '$');
if ~isempty(indr)
    name(indr) = 'D';
end
indr= find(name == '(');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == ')');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == '>');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == '<');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == ';');
if ~isempty(indr)
    name(indr) = '_';
end
if isempty(name)
    name = sprintf('noname%d',no);
    no = no + 1;
end
if (name(1)>='0') && (name(1)<='9')
    name = strcat('var_',name);
end

comment  = line(1:ind_deb(1));
value    = sscanf(line((ind_fin(end)+1):end),'%g');

% new format
if isempty(value)
    r = comment;
    nk = 100;
    while ~isempty(r) & nk > 0
        [s,r] = strtok(r,' ');
        v = sscanf(s,'%g');
        if ~isempty(v)
               value(end+1) = v;
        end
        nk = nk - 1;
    end  
    %value
end    



function  [name,value,comment,unit,no] = decode_line_var_alt(line,no);

name = line((line >= 58));
linef  = line(min(find((line >= 48) & (line <=57))):length(line));
if ~isempty(linef)
    value    = sscanf(linef(linef<58),'%g');
else
    value = [];
end
if isempty(value)
    value    = sscanf(line((max(find(line == '=')) + 1):end),'%g');
end
comment  = line((line >= 58) | (line <= 32));
if isempty(name) && ~isempty(value)
    name = sprintf('noname%d',no);
    no = no + 1;
elseif ~isempty(name) && ~isempty(value)
    name(name <= ' ') = '';
    indr= find(name == '=');
    if ~isempty(indr)
        name(indr) = '';
    end
    indr= find(name == '%');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == ',');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == '-');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == '+');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == '*');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == '.');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == ':');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == '/');
    if ~isempty(indr)
        name(indr) = 'o';
    end
    indr= find(name == '$');
    if ~isempty(indr)
        name(indr) = 'D';
    end
    indr= find(name == '(');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == ')');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == '>');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == '<');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == ';');
    if ~isempty(indr)
        name(indr) = '_';
    end
    indr= find(name == '~');
    if ~isempty(indr)
        name(indr) = '_';
    end
    if (name(1)>='0') && (name(1)<='9')
        name = strcat('var_',name);
    end
end
unit ='';
if length(name) > 63
    name = name(1:63);
end