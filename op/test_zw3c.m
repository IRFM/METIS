% test traduction auto
% url simple
[s,r,c]=zw3c('http://alpha.cad.cea.fr/index.html')
% method get + formulaire
url='http://alpha.cad.cea.fr/cgi-bin/logseq_cgi-bin.sh';
[s,r,c]=zw3c(url,'put','session','trait','nbligne','40')
% method post
url='http://tr.voila.fr/textrad';
%url='http://alpha.cad.cea.fr/cgi-bin/post-query';
[s,r,c]=zw3c(url,'post','template','mini','context','voila', ... 
             'lg','fr','template','Default','status','translate', ... 
             'direction','65544','kw','il fait beau et chaud');
	      
%value="524289">Anglais-Fran&ccedil;ais
%value="524292">Allemand-Fran&ccedil;ais
%value="262152">Fran&ccedil;ais-Allemand
%value="524290">Russe-Fran&ccedil;ais
%value="131080">Fran&ccedil;ais-Russe
%value="131073">Anglais-Russe
%value="65538">Russe-Anglais
%value="2097153">Anglais-Espagnol
%value="65540">Allemand-Anglais
%value="131076">Allemand-Russe
%value="262146">Russe-Allemand
%value="131088">Italien-Russe
