% appuie sur tous les boutons et capture les figures
function zpushtout(hf)
drawnow
try 
	hh = findobj(hf,'type','uicontrol','style','radio');
catch
   hh =[];
end
if isempty(hh)
    return
end

hmem =get(0,'children')


for k=1:length(hh)
  try
    cb= get(hh(k),'callback');
	 st = get(hh(k),'string');
	 if ~isempty(findstr(lower(st),'quit'))
	    ok = 0
	 elseif ~isempty(findstr(lower(st),'close'))
	    ok = 0
	 elseif ~isempty(findstr(lower(st),'ferm'))
	    ok = 0
	 elseif ~isempty(findstr(lower(st),'help'))
	    ok = 0
	 elseif ~isempty(findstr(lower(st),'aide'))
	    ok = 0
	 elseif ~isempty(findstr(lower(st),'canal'))
	    ok = 0
	 elseif ~isempty(findstr(lower(st),'annulation'))
	    ok = 0
	 elseif ~isempty(findstr(lower(st),'validation'))
	    ok = 0
	 elseif ~isempty(findstr(lower(st),'raz'))
	    ok = 0
	 else
	     ok = 1;
	 end
	 if ~isempty(cb) & (ok ==1)
	     set(hh(k),'value',1);
		  try
		    eval(cb)
			 pause(1)
			 hnew = get(0,'children');
			 hc   = setdiff(hnew,hmem);
			 if ~ isempty(hc)
			    for l =1:length(hc)
				     %figure(hc(l));
					  %drawnow
				     %zcapture(hc(l));
					  edit_ok =0;
					  tagc = get(hc(l),'tag');
					  if ~isempty(findstr(tagc,'editeur')) | ~isempty(findstr(tagc,'_mode')) | ...
					     ~isempty(findstr(tagc,'_profcmplx')) 
					      edit_ok =1;
					  end
					  if isempty(findobj(hc(l),'type','uicontrol','style','radio')) | edit_ok
				        figure(hc(l));
					     drawnow
				        %zcapture(hc(l));
					     delete(hc(l));
					  end 
				 end
			 end
		  end
	 end
  end
end
