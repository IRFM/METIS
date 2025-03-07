% extrait 1 temps d'une structure
function struct_out = zerod_get1t(struct_in,index)

if isstruct(struct_in)
	noms = fieldnames(struct_in);
	for k = 1:length(noms)
		struct_out.(noms{k}) = zerod_get1t(struct_in.(noms{k}) ,index);
	end
else
    sz = size(struct_in);
    if sz(1) > 1,
        struct_out = zeros([1,sz(2:end)]);
        struct_out(1,:) = struct_in(min(size(struct_in,1),index),:);
        if size(struct_in,1) < index
	    fprintf('X');
        end
    else
        struct_out = struct_in;
    end
end