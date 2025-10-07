%% Huy's Individual Code
% Assigned LOX_RP1 as Stage 2 Propellant
disp("running huy's code");

LOX_RP1.LOX_LCH4 = TwoStageLV([propellants.LOX_LCH4 propellants.LOX_RP1], DeltaV_Total, m_pl, delta, g);
LOX_RP1.LOX_LH2 = TwoStageLV([propellants.LOX_LH2 propellants.LOX_RP1], DeltaV_Total, m_pl, delta, g);
LOX_RP1.LOX_RP1 = TwoStageLV([propellants.LOX_RP1 propellants.LOX_RP1], DeltaV_Total, m_pl, delta, g);
LOX_RP1.Solid = TwoStageLV([propellants.Solid propellants.LOX_RP1], DeltaV_Total, m_pl, delta, g);
LOX_RP1.Storables = TwoStageLV([propellants.Storables propellants.LOX_RP1], DeltaV_Total, m_pl, delta, g);