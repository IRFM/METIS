% function to read URANIE/SYCOMORE tab 
function out = import_sycomore_tab(filename)

if nargin == 0
    [f,p] = uigetfile('*.dat');
    if isempty(f)
            return
    end
    filename = fullfile(p,f);
end


[data,delim,header] = importdata(filename,' ',4);
[void,out.name]     = strtok(data.textdata{1},' ');
out.name = strtrim(out.name);
[void,out.title]    = strtok(data.textdata{2},' ');
out.title = strtrim(out.title);
[void,out.date]     = strtok(data.textdata{3},' ');
out.date = strtrim(out.date);
[void,colname]      = strtok(data.textdata{4},' ');

k = 1;
while ~isempty(colname)
	colname = colname(2:end);
	[loc,colname] = strtok(colname,'|');
	if ~isempty(strfind(loc,'('))
		[fname,us]   = strtok(loc,'(');
		[unit,us]    = strtok(us,')');
		unit         = strtrim(unit(2:end));
		status       = strtrim(us(2:end));
	else
		fname =loc;
	end
	fname = formatname(strtrim(fname));
	out.(fname).value  = data.data(:,k);
	out.(fname).unit   = unit;
	out.(fname).status = status;
	
	k = k +1;
end

function name = formatname(name)

name(name <= ' ') = '';
indr= find(name == '=');
if ~isempty(indr)
    name(indr) = '_';
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
indr= find(name == '.');
if ~isempty(indr)
    name(indr) = '_';
end
indr= find(name == '/');
if ~isempty(indr)
    name(indr) = 'o';
end

