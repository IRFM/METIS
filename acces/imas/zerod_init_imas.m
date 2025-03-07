%==========================
%% READ DATA FROM IMAS UAL
%==========================
function z0dinput = zerod_init_imas(mode_exp,shot,gaz,temps,z0dinput)

prompt={'shot number :','Run :' ,'Occurrence:','Tokamak:','User:','Version:'};
def={'1','1','','','',''};
dlgTitle='Readig data from UAL';
lineNo=1;
answer=zinputdlg(prompt,dlgTitle,lineNo,def);
if isempty(answer)
  z0dinput = [];
  return
end
shot    = str2num(answer{1});
run     = str2num(answer{2});
occ     = answer{3};
ual_tok = answer{4};
ual_user= answer{5};
ual_ver = answer{6};

if ~isempty(ual_tok) && ~isempty(ual_tok(ual_tok > ' '))
  ual_occ = occ;
  clear occ
  occ.tokamak = ual_tok(ual_tok>' ');
  occ.user    = ual_user(ual_user> ' ');
  occ.dataversion = ual_ver(ual_ver> ' ');
  occ.occurrence = ual_occ(ual_occ> ' ');
end

[error_flag,z0dinput] = metis4imas(shot,run,occ,'metis_from_ual');
z0dinput.shot =  z0dinput.shot(1);
z0dinput.mode_exp = 117;
z0dinput.run = run;
if isstruct(occ)
  z0dinput.tokamak = occ.tokamak;
  z0dinput.user = occ.user;
  z0dinput.dataversion = occ.dataversion;
  z0dinput.occ = occ.occurrence;
else
  z0dinput.occ = occ;
end
z0dinput.shot = shot;
    

