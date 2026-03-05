%FUNTION TO READ the file consisting of ky values
function kys = read_ky_spectrum(simname)

	clear ky_values;
	fptr = fopen([simname, '/out.tglf.ky_spectrum']);
	frewind(fptr);

	% skip the first two lines
	fgetl(fptr);
	fgetl(fptr);
		

	% Initialize containers for the ky_vallues
	kys.ky_values = [];

	% Read file line by line
	tline = fgetl(fptr);

	while ischar(tline)
 	% Convert line to a numeric value and store it in the array
        	value = str2double(tline);
        	if ~isnan(value)
            		kys.ky_values(end+1) = value; %#ok<AGROW>
        	end
        	tline = fgetl(fptr);
    	end
    
    % Close the file
    fclose(fptr);

end






