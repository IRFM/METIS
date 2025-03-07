function zbpause

clk = clock;
%test du weekend
n   = datenum(clk(1),clk(2),clk(3),clk(4),clk(5),clk(6));
j   = lower(datestr(n,8));
switch j
case 'sat','sun'
   return
end

if (clk(4) > 7) & (clk(4) < 18)
    %pause_t = (18 - clk(4) -1) * 3600;
    %if pause_t > 0 
	 %		 pause(pause_t)
	 %end
	 [s,t]=unix('uptime');
	 if ((~isempty(t))&(s==0))
			indice=max(find(t==':'));
	 else
			indice=[];
	 end
	 if ~isempty(indice)
			charge=str2num(['[',t((indice+1):(length(t)-1)),']']);
			charge =charge(3);
	 else
			charge=inf;
	 end
	 
	 % calcul auto du seuil
	 [s,t] = unix('pset_info | wc -l');
	 if s ~= 0
	     seuil_cpu = 4
	 else
	     	try
	     		seuil_cpu = str2num(t) - 10;
			catch
	      	seuil_cpu = 4
	    	end
	 end
	 
	 fprintf('zbpause: %s  => %g / %d (%g)\n',datestr(now),charge,seuil_cpu,cputime);
	 while ((clk(4) > 7) & (clk(4) < 18) & (charge > seuil_cpu))
		 pause(300);
      clk = clock;
		[s,t]=unix('uptime');
		if ((~isempty(t))&(s==0))
			indice=max(find(t==':'));
		else
			indice=[];
		end
		if ~isempty(indice)
			charge=str2num(['[',t((indice+1):(length(t)-1)),']']);
			charge =charge(3);
		else
			charge=inf;
		end
	   fprintf('zbpause: %s  => %g / %d (%g)\n',datestr(now),charge,seuil_cpu,cputime);
	 end
end    
