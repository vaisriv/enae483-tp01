%% Huy's Individual Code
% Assigned LOX_RP1 as Stage 2 Propellant
disp("running huy's code");

%% LOX_LCH4
LOX_RP1.LOX_LCH4 = TwoStageLV([propellants.LOX_LCH4 propellants.LOX_RP1]);
LOX_RP1.LOX_LCH4 = LOX_RP1.LOX_LCH4.generate_trends(DeltaV_Total, m_pl, delta, g);
LOX_RP1.LOX_LCH4 = LOX_RP1.LOX_LCH4.save_mass_fig();
LOX_RP1.LOX_LCH4 = LOX_RP1.LOX_LCH4.save_cost_fig();

%% LOX_LH2
LOX_RP1.LOX_LH2 = TwoStageLV([propellants.LOX_LH2 propellants.LOX_RP1]);
LOX_RP1.LOX_LH2 = LOX_RP1.LOX_LH2.generate_trends(DeltaV_Total, m_pl, delta, g);
LOX_RP1.LOX_LH2 = LOX_RP1.LOX_LH2.save_mass_fig();
LOX_RP1.LOX_LH2 = LOX_RP1.LOX_LH2.save_cost_fig();

%% LOX_RP1
LOX_RP1.LOX_RP1 = TwoStageLV([propellants.LOX_RP1 propellants.LOX_RP1]);
LOX_RP1.LOX_RP1 = LOX_RP1.LOX_RP1.generate_trends(DeltaV_Total, m_pl, delta, g);
LOX_RP1.LOX_RP1 = LOX_RP1.LOX_RP1.save_mass_fig();
LOX_RP1.LOX_RP1 = LOX_RP1.LOX_RP1.save_cost_fig();

%% LOX_RP1
LOX_RP1.Solid = TwoStageLV([propellants.Solid propellants.LOX_RP1]);
LOX_RP1.Solid = LOX_RP1.Solid.generate_trends(DeltaV_Total, m_pl, delta, g);
LOX_RP1.Solid = LOX_RP1.Solid.save_mass_fig();
LOX_RP1.Solid = LOX_RP1.Solid.save_cost_fig();

%% Storables
LOX_RP1.Storables = TwoStageLV([propellants.Storables propellants.LOX_RP1]);
LOX_RP1.Storables = LOX_RP1.Storables.generate_trends(DeltaV_Total, m_pl, delta, g);
LOX_RP1.Storables = LOX_RP1.Storables.save_mass_fig();
LOX_RP1.Storables = LOX_RP1.Storables.save_cost_fig();
