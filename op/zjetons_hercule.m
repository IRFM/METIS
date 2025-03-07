% gestion des jetons nuit pour hercule
function je = zjetons_hercule

creation = 0;
je       = '';
try
    info = load('/usr/drfc/zineb/data/jetons_hercule.mat');
catch
    info ='';
end
if isempty(info)
	creation =1;
else
	clk =clock;
	if (clk(3) ~= info.clk(3)) & (clk(4) >= 8)
		% nouvelle journee
	   creation =1;
	end	
end

if creation == 1
	clk    = clock;
	jetons = {'23:00','00:00','01:00','02:00','03:00'};
	% mise en place du chien de garde
	% lancer par punch
 	%switch  datestr(datenum(clock),'ddd')
 	%case {'Tue','Wed','Thu'}
 	%      ! rsh rigel at 07:30 /usr/drfc/cgc/bin/zkiller
	%end      	
else
	clk    = info.clk;
	jetons = info.jetons;
end

if isempty(jetons)
	return
end

je     = jetons{1};
jetons = jetons(2:end);

save('/usr/drfc/zineb/data/jetons_hercule.mat','jetons','clk');
