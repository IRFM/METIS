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
	lp     = 'fr_en';
case 'allemand'
	lp     = 'fr_de';
otherwise
	lp     = 'fr_en';
end

% supression des cractere accentue
phrase = strrep(phrase,'�','e');
phrase = strrep(phrase,'�','E');
phrase = strrep(phrase,'�','e');
phrase = strrep(phrase,'�','e');
phrase = strrep(phrase,'�','e');
phrase = strrep(phrase,'�','c');
phrase = strrep(phrase,'�','a');
phrase = strrep(phrase,'�','a');
phrase = strrep(phrase,'�','a');
phrase = strrep(phrase,'�','u');
phrase = strrep(phrase,'�','u');
phrase = strrep(phrase,'�','u');
phrase = strrep(phrase,'�','i');

url = 'http://babel.altavista.com/tr';	
[s,r,c]=zw3c(url,'post','doit','done','tt','urltext','urltext',phrase,'lp',lp);
if s~=0
	rep =sprintf('Erreur de traduction (%d) :\n %s \n commande : \n %s \n',s,r,c);
	cr = s;		
else
	ind1  = findstr(r,'name="q">');
	ind2  = findstr(r,'</textarea>');
	if ~isempty(ind1)&~isempty(ind2)
		rep = r((ind1+9):(ind2-2));
	else
	   ind1  = findstr(r,'<td bgcolor=white>');
	   ind2  = findstr(r,'</td>');
		if ~isempty(ind1)&~isempty(ind2)
		   ind2 = min(ind2(ind2 >ind1)); 
		   rep = r((ind1+18):(ind2-2));
			rep(rep <32) =[];
	   else
		   rep =r;
		   cr = -999;
	   end
	end
end	      

ind = find(rep <32);
if ~isempty(ind)
	rep(ind) = ' ';
end

