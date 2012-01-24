% co2_co2
% This function will analyse the transcritical Rankine cycle for a
% power plant where the heat exchange fluid CO2 is run through a
% transcritical Rankine cycle to produce power
cycle_option = 'supercritical-loop';
%cycle_option = 'condensing-loop';
open_loop = 0;  % 0 for closed loop
loop_type = {'closed', 'open  '};
Trock = 273+235;    % Temperature of the hot rocks reservoir, oK
NP = 6;
figure(1); clf;
[T,s]=Ts_satcurve(200,0);
plot(s, T);
xlabel('Entropy, kJ/kg-^oK', 'FontSize', 12);
ylabel('Temperature, ^oK', 'FontSize', 12);
title(cycle_option);

hold on;
% Create the vectors to hold the properties at each point
v = zeros(6, 8);
%
% Point #1 - Assume the same state Pruess assumed for the starting point
switch cycle_option
    case 'condensing-loop'
        v(1,1) = 5.74; v(1,2) = 293; 
        % Find the complete state of the compressed liquid
        v1 = co2eqofstate('PTL', v(1,:));
    case 'supercritical-loop'
        v(1,1) = 8; v(1,2) = 320;
        v1 = co2eqofstate('PTG', v(1,:));
end

v(1,:)=v1;
% Mark the point on the chart
putpoint(v1(6), v1(2), '1');
%
% Point #2 - Bottom of Injection Hole
[v2, pinj, Tinj]=downhole(v1);
v(2,:)=v2;
figure(2); clf; plot(Tinj-273, -pinj, 'b'); axis([0 250 -70 0]); hold on;
ylabel('Pressure, MPa');
xlabel('Temperature, ^oC');
title('Pressure vs temperature down the injection hole')'
figure(1); putpoint(v(2,6),v(2,2),'2');
%
% Point #3 - Bottom of Production Hole
v(3,1)=v(2,1);
v(3,2)=Trock;
v3 = co2eqofstate('PTG', v(3,:));
v(3,:)=v3;
figure(1); putpoint(v(3,6),v(3,2),'3');
%
% Point #4 - Top of production hole
[v4, pprod, Tprod]=uphole(v(3,:));
v(4,:)=v4;
figure(2); plot(Tprod-273, -pprod, 'r')
% figure(1); putpoint(v(4,6),v(4,2),'4');



%
% Point %5 - After turbine
v(5,6) = v(4, 6);
v(5,1) = v(1, 1);
v(5,2) = v(1, 2);   % Trial value for temperature
%v(5,3) = co2prop('DV',v(5,2)); % Trial value for density
v5     = co2eqofstate('PTG', v(5,:));
v(5,3) = v5(3);     % Trial value for density
v5  = co2eqofstate('PS', v(5,:)); % Find the rest for given (p, s)
v(5, :)=v5;
figure(1); putpoint(v(5,6),v(5,2),'5');
%
% Point #6
switch cycle_option
    case 'condensing-loop'
        v(6,1) = v(5,1);
        v6  = co2eqofstate('PSATG', v(6,:));
    case 'supercritical-loop'
        v(6,:) = (v(5,:)+v(1,:))/2;
        v6  = v(6,:);
end
v(6,:) = v6;


% Redraw the cycle on the T-S diagram
figure(1); clf;
[T,s]=Ts_satcurve(200,0);
plot(s, T);
title(cycle_option);
xlabel('Entropy, kJ/kg-^oK');
ylabel('Temperature, ^oK');
hold on;
plot([v(6,6) v(1, 6)], [v(6,2) v(1, 2)], 'k--');
for i=1:5
    plot([v(i,6) v(i+1, 6)], [v(i,2) v(i+1, 2)], 'k--');
    putpoint(v(i,6),v(i,2),sprintf('%d',i));
end
% vtest = co2eqofstate('PS', v(4,:));
% figure(1); putpoint(vtest(6),vtest(2),'t');
% Calculate the cycle efficiency
qh = (v(3,6)-v(2,6))*(v(3,2)+v(2,2))/2;
qh2 = v(3,5)-v(2,5);
qh3 = v(4,5)-v(1,5);
qc = v(6,2)*(v(6,6)-v(1,6))+(v(5,6)-v(6,6))*(v(5,2)+v(6,2))/2;
qc2 = v(5,5)-v(1,5);
wt = v(4,5)-v(5,5); % Turbine work = h4-h5
if open_loop
    vin    = zeros(size(v(1,:)));
    vin(2) = 300;      % CO2 is available at this temperature as gas
    vin(3) = 1.7730;
    vin    = co2eqofstate('TD', vin);  % Properties for CO2 gas
    wp  = v(1,5)-vin(1,5);              % Pump work
    putpoint(vin(6), vin(2), '0');
else
    wp  = 0;
end
eff= (wt-wp)/qh*100;
eff2 = (qh2-qc2)/qh2;
ideal_eff=(1-v(1,2)/Trock)*100;
fprintf(1, '\n========================================================\n');
fprintf(1, '%s %s Simulation:\n', loop_type{open_loop+1}, cycle_option);
fprintf(1, 'Injection (top) conditions : %.1f MPa and %.0f ^oc\n', ...
    v(1,1), v(1,2)-273);
fprintf(1, 'Reservoir (downhole) conditions : %.1f MPa and %.0f ^oc\n', ...
    v(3,1), v(3,2)-273);
fprintf(1, 'Heat input into the cycle, qh(kJ/kg) = %.1f\n', qh);
fprintf(1, 'Turbine work, wt(kJ/kg)              = %.1f\n', wt);
fprintf(1, 'Pump work, wt(kJ/kg)                 = %.1f\n', wp);
fprintf(1, 'Thermal efficiency                   = %.1f%%\n', eff);
fprintf(1, 'Carnot efficiency                    = %.1f%%\n', ideal_eff);

text(-2.4, 500, sprintf('Thermal efficiency=%.1f%%', eff));



