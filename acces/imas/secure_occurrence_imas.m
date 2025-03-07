function occ_out = secure_occurrence_imas(occ)

if isempty(occ)
  occ_out = '';
elseif isnumeric(occ)
  occ = fix(occ);
  if le(occ,0)
    occ_out ='';
  else
    occ_out = num2str(occ);
  end
else
  occ = fix(str2num(occ));
  if le(occ,0)
    occ_out ='';
  else
    occ_out = num2str(occ);
  end
end
