%% Matthew's Individual Code
% Assigned LOX_LH2 as Stage 1 Propellant
disp("running matthew's code");

%% LOX_LCH4
LOX_LH2.LOX_LCH4 = TwoStageLV(propellants.LOX_LH2, propellants.LOX_LCH4);
LOX_LH2.LOX_LCH4 = LOX_LH2.LOX_LCH4.generate_trends(DeltaV_Total, m_pl, delta, g);
LOX_LH2.LOX_LCH4 = LOX_LH2.LOX_LCH4.save_mass_fig();
LOX_LH2.LOX_LCH4 = LOX_LH2.LOX_LCH4.save_cost_fig();

%% LOX_LH2
LOX_LH2.LOX_LH2 = TwoStageLV(propellants.LOX_LH2, propellants.LOX_LH2);
LOX_LH2.LOX_LH2 = LOX_LH2.LOX_LH2.generate_trends(DeltaV_Total, m_pl, delta, g);
LOX_LH2.LOX_LH2 = LOX_LH2.LOX_LH2.save_mass_fig();
LOX_LH2.LOX_LH2 = LOX_LH2.LOX_LH2.save_cost_fig();

%% LOX_RP1
LOX_LH2.LOX_RP1 = TwoStageLV(propellants.LOX_LH2, propellants.LOX_RP1);
LOX_LH2.LOX_RP1 = LOX_LH2.LOX_RP1.generate_trends(DeltaV_Total, m_pl, delta, g);
LOX_LH2.LOX_RP1 = LOX_LH2.LOX_RP1.save_mass_fig();
LOX_LH2.LOX_RP1 = LOX_LH2.LOX_RP1.save_cost_fig();

%% LOX_LH2
LOX_LH2.Solid = TwoStageLV(propellants.LOX_LH2, propellants.Solid);
LOX_LH2.Solid = LOX_LH2.Solid.generate_trends(DeltaV_Total, m_pl, delta, g);
LOX_LH2.Solid = LOX_LH2.Solid.save_mass_fig();
LOX_LH2.Solid = LOX_LH2.Solid.save_cost_fig();

%% Storables
LOX_LH2.Storables = TwoStageLV(propellants.LOX_LH2, propellants.Storables);
LOX_LH2.Storables = LOX_LH2.Storables.generate_trends(DeltaV_Total, m_pl, delta, g);
LOX_LH2.Storables = LOX_LH2.Storables.save_mass_fig();
LOX_LH2.Storables = LOX_LH2.Storables.save_cost_fig();
