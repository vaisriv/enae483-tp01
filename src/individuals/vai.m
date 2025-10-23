%% Vai's Individual Code
% Assigned Solid as Stage 2 Propellant
disp("running vai's code");

Solid.LOX_LCH4 = TwoStageLV([propellants.LOX_LCH4; propellants.Solid], DeltaV_Total, m_pl, payload_reqs, delta, g, d_lv);
Solid.LOX_LH2 = TwoStageLV([propellants.LOX_LH2; propellants.Solid], DeltaV_Total, m_pl, payload_reqs, delta, g, d_lv);
Solid.LOX_RP1 = TwoStageLV([propellants.LOX_RP1; propellants.Solid], DeltaV_Total, m_pl, payload_reqs, delta, g, d_lv);
Solid.Solid = TwoStageLV([propellants.Solid; propellants.Solid], DeltaV_Total, m_pl, payload_reqs, delta, g, d_lv);
Solid.Storables = TwoStageLV([propellants.Storables; propellants.Solid], DeltaV_Total, m_pl, payload_reqs, delta, g, d_lv);