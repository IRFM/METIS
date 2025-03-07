% ecalte un struture du workspace en sous structure
function zdeal(in)

nn = fieldnames(in);
for k =1:length(nn)
  out = getfield(in,nn{k});
  zassignin('base',nn{k},out);
end
