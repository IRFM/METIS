function [rep,cr] = ztraduit(phrase,langage)

cr =0;

if nargin <1
    rep ='';
    return
elseif isempty(phrase)
    rep ='';
    return
end
if nargin <2
   langage ='anglais';
end


switch lower(langage) 

case 'francais'
    rep =phrase;
    return
case 'anglais'
	value  = '65544';
	lp     = 'fr_en';
case 'allemand'
	value  = '262152';
	lp     = 'fr_de';
case 'russe'
   value  = '131080';
otherwise
   value  = '65544';
end


% supression des cractere accentue
phrase = strrep(phrase,'é','e');
phrase = strrep(phrase,'È','E');
phrase = strrep(phrase,'è','e');
phrase = strrep(phrase,'ë','e');
phrase = strrep(phrase,'ê','e');
phrase = strrep(phrase,'ç','c');
phrase = strrep(phrase,'à','a');
phrase = strrep(phrase,'â','a');
phrase = strrep(phrase,'ä','a');
phrase = strrep(phrase,'ù','u');
phrase = strrep(phrase,'û','u');
phrase = strrep(phrase,'ü','u');
phrase = strrep(phrase,'ï','i');

% method post
url='http://trans.voila.fr/textrad';
[s,r,c]=zw3c(url,'post','template','mini','context','voila', ... 
                 'lg','fr','template','Default','status','translate', ... 
                 'direction',value,'kw',phrase);

if s~=0
	rep =sprintf('Erreur de traduction (%d) :\n %s \n commande : \n %s \n', ...
	s,r,c);
	cr = s;		
else
	ind1  = findstr(r,'onFocus="clear()">');
	ind2  = findstr(r,'</textarea>');
	if ~isempty(ind1)&~isempty(ind2)
		rep = r((ind1+19):(ind2-3));
	else
		rep =r;
		cr = -999;
	end
end	      
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

ind = find(rep <32);
if ~isempty(ind)
	rep(ind) = ' ';
end
