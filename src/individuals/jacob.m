%% Jacob's Individual Code
% Assigned Storables as Stage 1 Propellant
disp("running jacob's code");

%% LOX_LCH4
Storables.LOX_LCH4 = TwoStageLV(propellants.Storables, propellants.LOX_LCH4);
Storables.LOX_LCH4 = Storables.LOX_LCH4.generate_trends(DeltaV_Total, m_pl, delta, g);
Storables.LOX_LCH4 = Storables.LOX_LCH4.save_mass_fig();
Storables.LOX_LCH4 = Storables.LOX_LCH4.save_cost_fig();

%% LOX_LH2
Storables.LOX_LH2 = TwoStageLV(propellants.Storables, propellants.LOX_LH2);
Storables.LOX_LH2 = Storables.LOX_LH2.generate_trends(DeltaV_Total, m_pl, delta, g);
Storables.LOX_LH2 = Storables.LOX_LH2.save_mass_fig();
Storables.LOX_LH2 = Storables.LOX_LH2.save_cost_fig();

%% LOX_RP1
Storables.LOX_RP1 = TwoStageLV(propellants.Storables, propellants.LOX_RP1);
Storables.LOX_RP1 = Storables.LOX_RP1.generate_trends(DeltaV_Total, m_pl, delta, g);
Storables.LOX_RP1 = Storables.LOX_RP1.save_mass_fig();
Storables.LOX_RP1 = Storables.LOX_RP1.save_cost_fig();

%% Storables
Storables.Solid = TwoStageLV(propellants.Storables, propellants.Solid);
Storables.Solid = Storables.Solid.generate_trends(DeltaV_Total, m_pl, delta, g);
Storables.Solid = Storables.Solid.save_mass_fig();
Storables.Solid = Storables.Solid.save_cost_fig();

%% Storables
Storables.Storables = TwoStageLV(propellants.Storables, propellants.Storables);
Storables.Storables = Storables.Storables.generate_trends(DeltaV_Total, m_pl, delta, g);
Storables.Storables = Storables.Storables.save_mass_fig();
Storables.Storables = Storables.Storables.save_cost_fig();
