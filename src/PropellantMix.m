% PropellantMix.m
classdef PropellantMix
        %PROPELLANTMIX engine propellant mixture container
        % class holding common propellant parameters

        properties
                % short name identifier (string)
                name string

                % display name for plot titles, legends, etc.
                displayname string

                % oxidizer to fuel mass ratio (ox/fuel)
                oxidizer_fuel_mix_mass_ratio double

                % specific impulse (s)
                specific_impulse double

                % thrust per motor for stage1 (N)
                thrust_per_motor_stage1 double

                % thrust per motor for stage2 (N)
                thrust_per_motor_stage2 double

                % engine exhaust diameter for stage1 (m)
                engine_exhaust_diameter_stage1 double

                % engine exhaust diameter for stage2 (m)
                engine_exhaust_diameter_stage2 double

                % chamber pressure for stage1 (pa)
                chamber_pressure_stage1 double

                % chamber pressure for stage2 (pa)
                chamber_pressure_stage2 double

                % nozzle expansion ratio for stage1 (ae/at)
                nozzle_expansion_ratio_stage1 double

                % nozzle expansion ratio for stage2 (ae/at)
                nozzle_expansion_ratio_stage2 double
        end

        methods
                function obj = PropellantMix(name, displayname, oxidizer_fuel_mix_mass_ratio, specific_impulse, thrust_per_motor_stage1, thrust_per_motor_stage2, engine_exhaust_diameter_stage1, engine_exhaust_diameter_stage2, chamber_pressure_stage1, chamber_pressure_stage2, nozzle_expansion_ratio_stage1, nozzle_expansion_ratio_stage2)
                        %constructor: assign all properties from inputs
                        %
                        % usage: p = PropellantMix(name, displayname, ox_f_ratio, isp, t1, t2, d1, d2, pc1, pc2, er1, er2)
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