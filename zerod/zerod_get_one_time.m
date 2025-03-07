% extrait 1 temps d'une structure
function struct_out = zerod_get_one_time(struct_in,time_in,time_out,method,extrap)

if isstruct(struct_in)
	noms = fieldnames(struct_in);
	for k = 1:length(noms)
		struct_out.(noms{k}) = zerod_get_one_time(struct_in.(noms{k}) ,time_in,time_out,method,extrap);
	end
else
    sz = size(struct_in);
    if (sz(1) > 1) && (sz(1) == length(time_in))
        struct_out = zeros([1,sz(2:end)]);
        for k=1:sz(2)
            struct_out(1,k) = interp1(time_in,struct_in(:,k),time_out(1),method,extrap);
        end
     else
        struct_out = struct_in;
    end
end