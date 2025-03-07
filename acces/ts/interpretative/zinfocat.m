function info = zinfocat(info,sdata,producteur,certin)

% cas signal simple
if nargin  == 4 
  nom = sdata;
  sdata =[];
  stu   = [];
  [st.cert,st.vers,st.date,st.heure,st.unix, ...
   st.uniy,st.uniz,st.nomdon,st.type]=tsbase_cert(certin);
  information =[];
  stu.commentaire = '';
  stu.information = st;
  sdata = setfield(sdata,nom,stu);
end

noms = fieldnames(sdata);
for k = 1:length(noms)
   nomc = noms{k};
   st = getfield(sdata,nomc);
   stc =[];
   if isfield(st,'information')
      stc = st.information;
   end
   if isfield(st,'commentaire')
      stc.com = st.commentaire;
   end
   if isfield(st,'data')
      stc.size = size(st.data);
   end
   stc.producteur = producteur;
   info =setfield(info,nomc,stc);
end
