%% Jacob's Individual Code
% Assigned Storables as Stage 2 Propellant
disp("running jacob's code");

Storables.LOX_LCH4 = TwoStageLV([propellants.LOX_LCH4; propellants.Storables], DeltaV_Total, m_pl, payload_reqs, delta, g, d_lv);
Storables.LOX_LH2 = TwoStageLV([propellants.LOX_LH2; propellants.Storables], DeltaV_Total, m_pl, payload_reqs, delta, g, d_lv);
Storables.LOX_RP1 = TwoStageLV([propellants.LOX_RP1; propellants.Storables], DeltaV_Total, m_pl, payload_reqs, delta, g, d_lv);
Storables.Solid = TwoStageLV([propellants.Solid; propellants.Storables], DeltaV_Total, m_pl, payload_reqs, delta, g, d_lv);
Storables.Storables = TwoStageLV([propellants.Storables; propellants.Storables], DeltaV_Total, m_pl, payload_reqs, delta, g, d_lv);