% script de test de zerodevolution
z0dinput = zerod_init(1,33850,2);
[zs,profil,z0dstruct] = zerodevolution([],z0dinput.option,z0dinput.cons.temps(1),z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
for k=2:length(z0dinput.cons.temps)
	[zs,profil,z0dstruct] = zerodevolution(z0dstruct,z0dinput.option,z0dinput.cons.temps(k),z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
	post.zerod = z0dstruct.zs;
	post.profil0d = z0dstruct.zs;
	post.z0dinput =z0dstruct.z0dinput;
	z0plotsc;
	drawnow
end