classdef PropellantMix
    %PROPELLANTMIX Engine Propellant Mixture

    properties
        oxidizer_fuel_mix_mass_ratio
        specific_impulse
        thrust_per_motor_stage1
        thrust_per_motor_stage2
        engine_exhaust_diameter_stage1
        engine_exhaust_diameter_stage2
        chamber_pressure_stage1
        chamber_pressure_stage2
        nozzle_expansion_ratio_stage1
        nozzle_expansion_ratio_stage2
    end

    methods
        function obj = PropellantMix(oxidizer_fuel_mix_mass_ratio, specific_impulse, thrust_per_motor_stage1, thrust_per_motor_stage2, engine_exhaust_diameter_stage1, engine_exhaust_diameter_stage2, chamber_pressure_stage1, chamber_pressure_stage2, nozzle_expansion_ratio_stage1, nozzle_expansion_ratio_stage2)
            %PROPELLANTMIX Construct an instance of this class
            obj.oxidizer_fuel_mix_mass_ratio = oxidizer_fuel_mix_mass_ratio;
            obj.specific_impulse = specific_impulse;
            obj.thrust_per_motor_stage1 = thrust_per_motor_stage1;
            obj.thrust_per_motor_stage2 = thrust_per_motor_stage2;
            obj.engine_exhaust_diameter_stage1 = engine_exhaust_diameter_stage1;
            obj.engine_exhaust_diameter_stage2 = engine_exhaust_diameter_stage2;
            obj.chamber_pressure_stage1 = chamber_pressure_stage1;
            obj.chamber_pressure_stage2 = chamber_pressure_stage2;
            obj.nozzle_expansion_ratio_stage1 = nozzle_expansion_ratio_stage1;
            obj.nozzle_expansion_ratio_stage2 = nozzle_expansion_ratio_stage2;
        end

        function [m, m_pr, m_o] = calculate_stage_masses(obj, DV, m_pl, delta, g)
            %CALCULATE_STAGE_MASSES Calculate stage masses
            r = exp(-DV/(g*obj.specific_impulse));
            lambda = r - delta;
            m_o = m_pl/lambda;
            m_in = delta*m_o;
            m_pr = m_o - m_pl - m_in;
            m = m_in + m_pr;
        end
    end
end