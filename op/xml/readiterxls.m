% first save as txt, then as .xls

filename='cond.xls';
[a,b]=xlsread(filename);

rz=find(strcmp(b(:,2),'R, m'))+1  ; 

innerr = a(rz(1)+[0:min(find(isnan(a(rz(1):end,2)))-2)],2);
innerz = a(rz(1)+[0:min(find(isnan(a(rz(1):end,3)))-2)],3);

outerr = a(rz(3):end,2);
outerz = a(rz(3):end,3);

supportr = a(rz(2)+[0:min(find(isnan(a(rz(2):end,2)))-2)],2);
supportz = a(rz(2)+[0:min(find(isnan(a(rz(2):end,3)))-2)],3);

filename = 'pf.xls';
[a,b]=xlsread(filename);

indpf1 = [8:19]-2;
indpf=indpf1;
coils = b(indpf+2,2);
rcoils = a(indpf,3);
zcoils = a(indpf,4);
drcoils = a(indpf,5);
dzcoils = a(indpf,6);
ntcoils = a(indpf,7);

indpf2 = indpf1+16;
indpf=indpf2;
imaxsod=a(indpf,4)*1e3;
imaxeob=a(indpf,5)*1e3;
bmaxsod=a(indpf,6);
bmaxeob=a(indpf,7);

indpf3 = indpf2 + 32;
indpf=indpf3;
umax = a(indpf3,3)*1e3;

pstepmax = 60e6;
ppfmax = 100e6;

filename='fw.xls';
[a,b]=xlsread(filename);

rlim = a(6:end,2);
zlim = a(6:end,3);

