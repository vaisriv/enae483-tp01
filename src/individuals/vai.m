%% Vai's Individual Code
% Assigned Solid as Stage 1 Propellant
disp("running vai's code");

%% LOX_LCH4
Solid.LOX_LCH4 = TwoStageLV(propellants.Solid, propellants.LOX_LCH4);
Solid.LOX_LCH4 = Solid.LOX_LCH4.generate_trends(DeltaV_Total, m_pl, delta, g);
Solid.LOX_LCH4 = Solid.LOX_LCH4.save_mass_fig();
Solid.LOX_LCH4 = Solid.LOX_LCH4.save_cost_fig();

%% LOX_LH2
Solid.LOX_LH2 = TwoStageLV(propellants.Solid, propellants.LOX_LH2);
Solid.LOX_LH2 = Solid.LOX_LH2.generate_trends(DeltaV_Total, m_pl, delta, g);
Solid.LOX_LH2 = Solid.LOX_LH2.save_mass_fig();
Solid.LOX_LH2 = Solid.LOX_LH2.save_cost_fig();

%% LOX_RP1
Solid.LOX_RP1 = TwoStageLV(propellants.Solid, propellants.LOX_RP1);
Solid.LOX_RP1 = Solid.LOX_RP1.generate_trends(DeltaV_Total, m_pl, delta, g);
Solid.LOX_RP1 = Solid.LOX_RP1.save_mass_fig();
Solid.LOX_RP1 = Solid.LOX_RP1.save_cost_fig();

%% Solid
Solid.Solid = TwoStageLV(propellants.Solid, propellants.Solid);
Solid.Solid = Solid.Solid.generate_trends(DeltaV_Total, m_pl, delta, g);
Solid.Solid = Solid.Solid.save_mass_fig();
Solid.Solid = Solid.Solid.save_cost_fig();

%% Storables
Solid.Storables = TwoStageLV(propellants.Solid, propellants.Storables);
Solid.Storables = Solid.Storables.generate_trends(DeltaV_Total, m_pl, delta, g);
Solid.Storables = Solid.Storables.save_mass_fig();
Solid.Storables = Solid.Storables.save_cost_fig();
