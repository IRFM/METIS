function rapfile=z0rapport(rapfile)

if nargin ==0
	info = evalin('base','post.z0dinput');
	rapfile = sprintf('z0d_%s_%d',info.machine,fix(info.shot));
end
warning off
try
	delete(sprintf('%s.ps',rapfile));
	delete(sprintf('%s.pdf',rapfile));
end

% demande du temps initial et final
% choix de la langue
langue =  lower(getappdata(0,'langue_cronos'));
temps  = evalin('base','post.zerod.temps');
tdeb = num2str(min(temps));
tfin = num2str(max(temps));

switch langue
case 'francais'
      prompt={'temps de debut (s):','temps de fin (s):'};
      def={tdeb,tfin};
      dlgTitle='Base temps pour le 0D';
otherwise
      prompt={'begin time (s):','end time (s):'};
      def={tdeb,tfin};
      dlgTitle='0D time slices';
end
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
if isempty(answer)
   	return
end
tdeb  = str2num(answer{1});
tfin  = str2num(answer{2});


hfl = get(0,'children');
zuicreefunform('z0dall','z0dinput.option',1,1);
fp1(rapfile,1);

evalin('base','z0plotsc');
tlim(gcf,tdeb,tfin);
fp1(rapfile);

evalin('base','z0plotp');
tlim(gcf,tdeb,tfin);
fp1(rapfile);

evalin('base','z0plote');
tlim(gcf,tdeb,tfin);
fp1(rapfile);

evalin('base','z0plotj');
tlim(gcf,tdeb,tfin);
fp1(rapfile);

evalin('base','z0plotn');
tlim(gcf,tdeb,tfin);
fp1(rapfile);

evalin('base','z0plott');
tlim(gcf,tdeb,tfin);
fp1(rapfile);

evalin('base','z0plotc');
tlim(gcf,tdeb,tfin);
fp1(rapfile);

evalin('base','z0ploteq');
tlim(gcf,tdeb,tfin);
fp1(rapfile);

evalin('base','z0plotlh');
tlim(gcf,tdeb,tfin);
fp1(rapfile);

evalin('base','z0plotgeo');
tlim(gcf,tdeb,tfin);
fp1(rapfile);

evalin('base','z0plotconv');
tlim(gcf,tdeb,tfin);
fp1(rapfile);

% conversion
[s,t]=unix(sprintf('ps2pdf %s.ps %s.pdf',rapfile,rapfile));
if s ~=0
	error(t);
else
	delete(sprintf('%s.ps',rapfile))
end
% fin
warning on
hfln = get(0,'children');
for k=1:length(hfln)
    if ishandle(hfln(k)) && ~any(hfl == hfln(k))
            close(hfln(k));
    end
end



function fp1(nom,first)


hf = gcf;
figure(hf);
pos = get(hf,'position');

% selon la version de matlab et de l'OS
try
   if nargin < 2
	set(hf,'PaperUnits','centimeters','PaperOrientation','portrait','PaperPosition',[1.5 1.5 18 23], ...
        	'PaperPositionMode','manual','PaperSize',[20.984 29.6774],'PaperType','A4');
   else
	set(hf,'PaperUnits','centimeters','PaperOrientation','portrait', ...
        	'PaperPositionMode','auto','PaperType','A4');
   end
catch
   if nargin < 2
	set(hf,'PaperUnits','centimeters','PaperOrientation','portrait','PaperPosition',[1.5 1.5 18 23], ...
        	'PaperPositionMode','manual','PaperSize',[20.984 29.6774],'PaperType','A4');
   else
	set(hf,'PaperUnits','centimeters','PaperOrientation','portrait', ...
        	'PaperPositionMode','auto','PaperType','A4');
   end
end


drawnow;
legend;
drawnow
if nargin < 2
	print(hf,'-append','-dpsc2','-noui',nom);
	%print('-append','-dpsc2','-noui',sprintf('-f%d',hf),nom);
	%print('-append','-dpsc2',sprintf('-f%d',hf),nom);
else
	print(hf,'-append','-dpsc2',nom);
	%print('-append','-dpsc2',sprintf('-f%d',hf),nom);
end
%ls('-l',strcat(nom,'*'))
%close(hf);



function tlim(hf,tdeb,tfin)

ha = findobj(hf,'type','axes');
set(ha,'xlim',[tdeb,tfin],'ylimmode','auto','yscale','linear');
drawnow
for k=1:length(ha)
	hc = ha(k);
	ylim = get(hc,'ylim');
	if ylim(1) > 0
		ylim(1) =0;
	end
	set(hc,'ylim',ylim);
end





