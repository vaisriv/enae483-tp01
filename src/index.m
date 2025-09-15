DV_Total = 12.3;
m_pl = 26000;

propellants = struct(...
    LOX_LCH4 = PropellantMix(3.6, 327, 2.26, 0.745, 2.4, 1.5, 35.16, 10.1, 34.34, 45), ...
    LOX_LH2 = PropellantMix(6.03, 366, 1.86, 0.099, 2.4, 2.15, 20.64, 4.2, 78, 84), ...
    LOX_RP1 = PropellantMix(2.72, 311, 1.92, 0.061, 3.7, 0.92, 25.8, 6.77, 37, 14.5), ...
    Solid = PropellantMix(1, 269, 4.5, 2.94, 6.6, 2.34, 10.5, 5, 16, 56), ...
    Storable = PropellantMix(2.67, 285, 1.75, 0.067, 1.5, 1.13, 15.7, 14.7, 26.2, 81.3) ...
);

delta1 = 0.08;
delta2 = 0.08;

g = 9.82;

X = linspace(0, 1, 100);
Y1 = zeros(1, 100);
Y2 = zeros(1, 100);
for i=1:100
    Y1(i) = propellants.LOX_LCH4.calculate_stage_masses(DV_Total*X(i), m_pl, delta1, g);
    Y2(i) = propellants.Solid.calculate_stage_masses(DV_Total*(1-X(i)), m_pl, delta2, g);
end

figure;
hold on;
plot(DV_Total*X, Y1, Color="g");
plot(DV_Total*X, Y2, Color="b");
hold off;