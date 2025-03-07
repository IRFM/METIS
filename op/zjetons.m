% gestion des jetons nuit pour rigel
function je = zjetons

creation = 0;
je       = '';
try
    info = load('/usr/drfc/cgc/cgc_data/zineb/data/jetons.mat');
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
%	jetons = {'22:00','22:30','23:00','23:30', ...
%	          '00:00','00:30','01:30','02:00', ...
%	          '02:30','03:00','03:30','04:00'};
	jetons = {'22:00','23:00'};
	% mise en place du chien de garde
	% lancer par punch
% 	switch  datestr(datenum(clock),'ddd')
% 	case {'Mon','Tue','Wed','Thu'}
% 	      ! rsh rigel at 07:30 /usr/drfc/cgc/bin/zkiller
%    end      	
else
	clk    = info.clk;
	jetons = info.jetons;
end

if isempty(jetons)
	return
end

je     = jetons{1};
jetons = jetons(2:end);

save('/usr/drfc/cgc/cgc_data/zineb/data/jetons.mat','jetons','clk');
