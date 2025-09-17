classdef PropellantMix
        %PROPELLANTMIX Engine Propellant Mixture

        properties
                name string
                displayname string
                oxidizer_fuel_mix_mass_ratio double
                specific_impulse double
                thrust_per_motor_stage1 double
                thrust_per_motor_stage2 double
                engine_exhaust_diameter_stage1 double
                engine_exhaust_diameter_stage2 double
                chamber_pressure_stage1 double
                chamber_pressure_stage2 double
                nozzle_expansion_ratio_stage1 double
                nozzle_expansion_ratio_stage2 double
        end

        methods
                function obj = PropellantMix(name, displayname, oxidizer_fuel_mix_mass_ratio, specific_impulse, thrust_per_motor_stage1, thrust_per_motor_stage2, engine_exhaust_diameter_stage1, engine_exhaust_diameter_stage2, chamber_pressure_stage1, chamber_pressure_stage2, nozzle_expansion_ratio_stage1, nozzle_expansion_ratio_stage2)
                        %PROPELLANTMIX Construct an instance of this class
                        obj.name = name;
                        obj.displayname = displayname;
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
        end
end