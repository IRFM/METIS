%testing new zsigmavfusion2 function.



% comparing in linear space

tikeV1 = linspace(1, 100e3);
[dd_p1,dd_n1,dt1,dhe31, ~, ~, ~]=zsigmavfusion(tikeV1);

tikeV2 = linspace(1, 1000e3);
[dd_p2,dd_n2,dt2,dhe32, ~, ~, ~]=zsigmavfusion2(tikeV2);

figure (1)
clf;
plot(tikeV1, dd_p1.sv, '-b', 'LineWidth', 2); hold on;
plot(tikeV1, dd_n1.sv, '-g', 'LineWidth', 2);
plot(tikeV1, dt1.sv, '-r', 'LineWidth', 2);
plot(tikeV1, dhe31.sv, '-c', 'LineWidth', 2);

plot(tikeV2, dd_p2.sv, 'ob', 'LineWidth', 2);
plot(tikeV2, dd_n2.sv, 'og', 'LineWidth', 2);
plot(tikeV2, dt2.sv, 'or', 'LineWidth', 2);
plot(tikeV2, dhe32.sv, 'oc', 'LineWidth', 2);

legend('dd\_p1', 'dd\_n1', 'dt1', 'dhe31', 'dd\_p2', 'dd\_n2', 'dt2', 'dhe32', 'Location', 'best');
xlabel('Temperature (eV)');
ylabel('Reaction Rate (cm^2/s)');
title('Fusion Cross Section vs Energy');
grid on;



% comparing in log space

tikeV1 = logspace(1, 5);
[dd_p1,dd_n1,dt1,dhe31, ~, ~, ~]=zsigmavfusion(tikeV1);

tikeV2 = logspace(1, 6);
[dd_p2,dd_n2,dt2,dhe32, ~, ~, ~]=zsigmavfusion2(tikeV2);

figure (2)
clf;
loglog(tikeV1, dd_p1.sv, '-b', 'LineWidth', 2); hold on;
loglog(tikeV1, dd_n1.sv, '-g', 'LineWidth', 2);
loglog(tikeV1, dt1.sv, '-r', 'LineWidth', 2);
loglog(tikeV1, dhe31.sv, '-c', 'LineWidth', 2);

loglog(tikeV2, dd_p2.sv, 'ob', 'LineWidth', 2);hold on;
loglog(tikeV2, dd_n2.sv, 'og', 'LineWidth', 2);
loglog(tikeV2, dt2.sv, 'or', 'LineWidth', 2);
loglog(tikeV2, dhe32.sv, 'oc', 'LineWidth', 2);

legend('dd\_p1', 'dd\_n1', 'dt1', 'dhe31', 'dd\_p2', 'dd\_n2', 'dt2', 'dhe32', 'Location', 'best');
xlabel('Temperature (eV)');
ylabel('Reaction Rate (cm^2/s)');
title('Reaction Rate vs Temperature');
grid on;