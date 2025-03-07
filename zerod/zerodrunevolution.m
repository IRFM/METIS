% script pour le calcul de metis avec zerodevolution
% securite
if isfield(z0dinput.exp0d,'Rsepa') && isempty(z0dinput.exp0d.Rsepa)
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Rsepa');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'Zsepa');
end
if isfield(z0dinput.exp0d,'XDURx') && isempty(z0dinput.exp0d.XDURx)
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURx');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURt');
	z0dinput.exp0d = rmfield(z0dinput.exp0d,'XDURv');
end
% initialisation
cons1t  = zerod_get1t(z0dinput.cons,1);
geo1t   = zerod_get1t(z0dinput.geo,1);
exp0d1t = zerod_get1t(z0dinput.exp0d,1);
[zs,profil,z0dstruct] = zerodevolution([],z0dinput.option,z0dinput.cons.temps(1),cons1t,geo1t,z0dinput.exp0d,exp0d1t,exp0d1t);
hdlg = z0dpatience('evolution');
% boucle sur les temps :
for k_index = 2:length(z0dinput.cons.temps)
        z0dpatience(k_index/(length(z0dinput.cons.temps)+1));
	cons1t  = zerod_get1t(z0dinput.cons,k_index);
	geo1t   = zerod_get1t(z0dinput.geo,k_index);
	exp0d1t = zerod_get1t(z0dinput.exp0d,k_index);
	[zs,profil,z0dstruct] = zerodevolution(z0dstruct,z0dinput.option,z0dinput.cons.temps(k_index),cons1t,geo1t,z0dinput.exp0d,exp0d1t,exp0d1t);
end
z0dpatience(1);

post.z0dinput = z0dstruct.z0dinput;
post.zerod    = z0dstruct.zs;
post.profil0d = z0dstruct.profil;
delete(hdlg)