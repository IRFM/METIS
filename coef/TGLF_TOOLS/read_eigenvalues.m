% To be generalized for any number of eigenmodes
function eigenvalues = read_eigenvalues(simname)

    clear eigenvalues;
    fptr = fopen([simname, '/out.tglf.eigenvalue_spectrum']);
    frewind(fptr);

    % skip the first two lines
    fgetl(fptr);
    fgetl(fptr);

    % Initialize containers for the eigenvalues
    eigenvalues.gamma1 = [];
    eigenvalues.freq1 = [];
    eigenvalues.gamma2 = [];
    eigenvalues.freq2 = [];

    % Read file line by line
    tline = fgetl(fptr);

    while ischar(tline)
        % Parse numerical values from the line
        nums = sscanf(tline, '%f');
        
        % Check the number of elements and store accordingly
        if numel(nums) == 4
            eigenvalues.gamma1(end+1) = nums(1);
            eigenvalues.freq1(end+1) = nums(2);
            eigenvalues.gamma2(end+1) = nums(3);            
            eigenvalues.freq2(end+1) = nums(4);
        end
        
        % Read the next line
        tline = fgetl(fptr);
    end
    
    % Close the file
    fclose(fptr);

end
