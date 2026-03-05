function outrun = read_tglf_outrun_3sp(simname)

% This function reads fluxes from the TGLF output namelist file out.tglf.run as provided by the "simname" input string (given with full path and name) into the matlab structure "outrun"

%simname = '/Home/RK278477/codes/metis/coef/TGLF_TOOLS/runs/tglf_tmp_radial2';
clear outrun;
outrun = struct();
fptr = fopen([simname, '/out.tglf.run']);

frewind(fptr);

for ij = 1:100
	aaa = fgetl(fptr);
	if length(aaa) > 4
		if strcmp(aaa(1:4), 'elec')
			outrun.elec_flx = str2num(aaa(6:end));
		end

		if strcmp(aaa(1:4), 'ion1')
			outrun.ion1_flx = str2num(aaa(6:end));
		end

		if strcmp(aaa(1:4), 'ion2')
			outrun.ion2_flx = str2num(aaa(6:end));
		end
	end
end
fclose(fptr);


%if ~isfield(outrun, 'ion2_flx')
%	outrun.ion2_flx = outrun.ion1_flx*0.0;
%end

outrun.Ge_gB = outrun.elec_flx(1);
outrun.G1_gB = outrun.ion1_flx(1);
outrun.G2_gB = outrun.ion2_flx(1);
outrun.Qe_gB = outrun.elec_flx(2);
outrun.Q1_gB = outrun.ion1_flx(2);
outrun.Q2_gB = outrun.ion2_flx(2);

end

