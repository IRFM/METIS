function data = read_sycomore_output_file(filename)

data = struct;
fid = fopen(filename);
if fid < 0
    fprintf('uanble to open file %s\n',filename);
    return
end

% boucle de lecture
tk = 1;
no = 1;
for k=1:10000
    line = fgetl(fid);
    if ~ischar(line)
      break
    end
    inde = find(line  == '=');
    if ~isempty(inde) & isempty(strfind(line,'=='))
	if length(inde) < 3
	    inde = max(inde);
	    name = line(1:(inde - 1));
	    info = line((inde + 1):end);
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
	    
	    num = sscanf(info,'%g');
	    if ~isempty(num)
		  data.(name) = num;
	    else
		  data.(name) = info;
	  
	    end
	 else
	    data.(sprintf('title%d',tk)) = line;
	    tk = tk + 1;
	 end
    elseif strfind(line,'---------------')
	  [name,s] = strtok(line,'---------------');	    
	  num = sscanf(s((length('---------------')+1):end),'%g');
	  if ~isempty(num)
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
	      if isempty(name)
		    data.(sprintf('nomane%d',no)) = num;
		    no = no + 1;
	      else
		    data.(name) = num;
	      end
	  end
    else
	    data.(sprintf('title%d',tk)) = line;
	    tk = tk + 1;
    end
end

fclose(fid);
