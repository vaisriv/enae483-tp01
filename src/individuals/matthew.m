%% Matthew's Individual Code
% Assigned LOX_LH2 as Stage 2 Propellant
disp("running matthew's code");

LOX_LH2.LOX_LCH4 = TwoStageLV([propellants.LOX_LCH4; propellants.LOX_LH2], DeltaV_Total, m_pl, delta, g, d_lv);
LOX_LH2.LOX_LH2 = TwoStageLV([propellants.LOX_LH2; propellants.LOX_LH2], DeltaV_Total, m_pl, delta, g, d_lv);
LOX_LH2.LOX_RP1 = TwoStageLV([propellants.LOX_RP1; propellants.LOX_LH2], DeltaV_Total, m_pl, delta, g, d_lv);
LOX_LH2.Solid = TwoStageLV([propellants.Solid; propellants.LOX_LH2], DeltaV_Total, m_pl, delta, g, d_lv);
LOX_LH2.Storables = TwoStageLV([propellants.Storables; propellants.LOX_LH2], DeltaV_Total, m_pl, delta, g, d_lv);