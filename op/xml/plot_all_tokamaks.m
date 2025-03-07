%% plot all tokamaks

clear all
close all

tokamaks = {'tcv','iter','aug','mast','globus'}

for i=1:length(tokamaks)
 figure(i)
 plot_tokamak(tokamaks{i})
 title(['DINA-CH General Tokamak Defintion - ' upper(tokamaks{i})])
 drawnow
end

return

for i=1:length(tokamaks)
 figure(i)
 print -Pjo
end
