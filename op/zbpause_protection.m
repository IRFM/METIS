function zbpause_protection

clk = clock;
%test du weekend
n   = datenum(clk(1),clk(2),clk(3),clk(4),clk(5),clk(6));
j   = lower(datestr(n,8));
switch j
case 'sat','sun'
   return
end

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

% anti pompage
seuil_cpu = 1.1 .* seuil_cpu ;
	 
fprintf('zbpause: %s  => %g / %d (%g)\n',datestr(now),charge,seuil_cpu,cputime);
while (charge > seuil_cpu)
		 pause(30);
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
			charge = 0;
		end
	   fprintf('zbpause: %s  => %g / %d (%g)\n',datestr(now),charge,seuil_cpu,cputime);
end    
