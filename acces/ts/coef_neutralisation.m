function [coef_h,coef_d,energy] = coef_neutralisation(idisplay)
%% ---------------------------------------------------------------
%% PURPOSE: CALCULATE THE NEUTRALIZATION EFFICIENCY AS A FUNCTION
%% OF THE BEAM ENERGY
%% ---------------------------------------------------------------
%% OPTIONAL ARGUMENT: IDISPLAY = 1 TO PLOT THE NEUTRALIZATION
%%                    COEFFICIENT VERSUS THE BEAM ENERGY
%% ---------------------------------------------------------------
%% REFERENCE = "LA FUSION THERMONUCLEAIRE CONTROLEE PAR
%% CONFINEMENT MAGNETIQUE", COLLECTION CEA, P167, FIG. 3.2.6
%% ---------------------------------------------------------------

if nargin < 1
  idisplay = 0;
end

%% ENERGY RANGE BETWEEN 0 AND 250 KEV
energy = linspace(0,250,26);

%% NEUTRALIZATION COEFFICIENT FOR HYDROGEN
coef_h(1)  = 11.4;
coef_h(2)  = 10.8;
coef_h(3)  = 10.1;
coef_h(4)  = 8.9;
coef_h(5)  = 7.6;
coef_h(6)  = 6.4;
coef_h(7)  = 5.3;
coef_h(8)  = 4.4;
coef_h(9)  = 3.65;
coef_h(10) = 3;
coef_h(11) = 2.5;
coef_h(12) = 2.05;
coef_h(13) = 1.75;
coef_h(14) = 1.35;
coef_h(15) = 1.15;
coef_h(16) = 1.;
coef_h(17) = 0.8;
coef_h(18) = 0.65;
coef_h(19) = 0.5;
coef_h(20) = 0.4;
coef_h(21) = 0.3;
coef_h(22) = 0.2;
coef_h(23) = 0.17;
coef_h(24) = 0.15;
coef_h(25) = 0.10;
coef_h(26) = 0.10;

coef_h = coef_h * 0.9/11.15;

%% NEUTRALIZATION COEFFICIENT FOR DEUTERIUM
coef_d = coef_h * sqrt(2.);

%% DISPLAY NEUTRALIZATION COEFFICIENT VERSUS BEAM ENERGY
if idisplay == 1
  
  figure
  
  linewidth = 2;
  fontsize  = 15;

  clf
  h=axes;
  set(h,'FontSize',fontsize)

  plot(energy,coef_h,'LineWidth',linewidth)
  hold on
  plot(energy,coef_d,'r--','LineWidth',linewidth)
  grid on
  xlabel('Beam energy (keV)','FontSize',fontsize)
  ylabel('Neutralization coefficient (-)','FontSize',fontsize)
  h=legend('Hydrogen','Deuterium','Location','NorthEast');
  set(h,'FontSize',fontsize)

end

