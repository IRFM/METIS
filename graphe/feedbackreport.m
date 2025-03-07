% display some useful figures for DIIID feedback simulations
% syntax : simply type feedbackreport
% load the input file in advance to display the feed forward NBI power


figure
if exist('jeux1')
   plot(data.gene.temps,data.gene.paddidn,'b-.',data.gene.temps,sum(jeux1.data.cons.idn')','g-')
   legend('result','feed forward')
else
   plot(data.gene.temps,data.gene.paddidn,'b-.')
end
xlabel('Time (s)')
ylabel('NBI power (W)')

figure
plot(data.gene.temps,data.gene.qmin,'b-.',data.gene.temps,data.cons.asser.qmin,'g-')
xlabel('Time (s)')
ylabel('qmin')
legend('result','target')

