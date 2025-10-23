%% Jack's Individual Code
% Assigned LOX_LCH4 as Stage 2 Propellant
disp("running jack's code");

LOX_LCH4.LOX_LCH4 = TwoStageLV([propellants.LOX_LCH4; propellants.LOX_LCH4], DeltaV_Total, m_pl, delta, g, d_lv);
LOX_LCH4.LOX_LH2 = TwoStageLV([propellants.LOX_LH2; propellants.LOX_LCH4], DeltaV_Total, m_pl, delta, g, d_lv);
LOX_LCH4.LOX_RP1 = TwoStageLV([propellants.LOX_RP1; propellants.LOX_LCH4], DeltaV_Total, m_pl, delta, g, d_lv);
LOX_LCH4.Solid = TwoStageLV([propellants.Solid; propellants.LOX_LCH4], DeltaV_Total, m_pl, delta, g, d_lv);
LOX_LCH4.Storables = TwoStageLV([propellants.Storables; propellants.LOX_LCH4], DeltaV_Total, m_pl, delta, g, d_lv);