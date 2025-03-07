% ok = 1 si TS n'est pas en exploitation
function ok = zetatTS

ok   = 0;
clk  = clock;
jour = sprintf('%4.4d%2.2d%2.2d',clk(1),clk(2),clk(3));
requete_nb = strrep(['SELECT ifnull(min (date_journee),0) from exploitation ', ...
	           'WHERE date_journee > $JOUR$ AND type_journee=''Exploitation'''],'$JOUR$',jour);

requete_1 = strrep(['SELECT ifnull(min (date_journee),0) from exploitation ', ...
	           'WHERE date_journee = $JOUR$ AND type_journee=''Exploitation'''],'$JOUR$',jour);
				  
[cr,temps,requete,texte,txt,nb_para,nom_para,rep]=matlsql('','',requete_1);
if rep ~= 0
    return
end

[cr,temps,requete,texte,txt,nb_para,nom_para,rep]=matlsql('','',requete_nb);
if rep == 0 
  ok = 1
  return
end

aa = fix(rep /10000);
rep = rep - aa * 10000;
mm = fix(rep /100);
jj = fix(rep - mm * 100);

nb   = (aa -clk(1)) * 365.25 + (mm - clk(2)) * 30.5 + (jj - clk(3));

if nb > 1
  ok = 1;
  fprintf('Reprise de l''exploitation de TS dans %d jours (+/- 1)\n',ceil(nb));
end
