% decompact une donnees cronos
function out = zuncompact(in)

if isstruct(in)
	noms = fieldnames(in);
	for k=1:length(noms)
		out.(noms{k}) =  zuncompact(in.(noms{k}));
	end
elseif issparse(in)
	out = full(in);
elseif isstr(in)
	if ~isempty(findstr(in,'ones(')) || ~isempty(findstr(in,'zeros('))
		warning off
		out = eval(in,'in');
		warning on
	else
		out = in;
	end
else
	out = in;
end