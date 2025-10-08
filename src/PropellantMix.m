% PropellantMix.m
classdef PropellantMix
        %PROPELLANTMIX engine propellant mixture container
        % class holding common propellant parameters

        properties
                % short name identifier (string)
                name string = "Blank"

                % display name for plot titles, legends, etc.
                displayname string = "Blank"

                % oxidizer mass ratio (ox/fuel) [-]
                oxidizer_mix_mass_ratio double = 0

                % fuel mass ratio (ox/fuel) [-]
                fuel_mix_mass_ratio double = 0

                % specific impulse [s]
                specific_impulse double = 0

                % thrust per motor for stage1 [N]
                thrust_per_motor_stage1 double = 0

                % thrust per motor for stage2 [N]
                thrust_per_motor_stage2 double = 0

                % engine exhaust diameter for stage1 [m]
                engine_exhaust_diameter_stage1 double = 0

                % engine exhaust diameter for stage2 [m]
                engine_exhaust_diameter_stage2 double = 0

                % chamber pressure for stage1 [Pa]
                chamber_pressure_stage1 double = 0

                % chamber pressure for stage2 [Pa]
                chamber_pressure_stage2 double = 0

                % nozzle expansion ratio for stage1 (A_e/A_t) [-]
                nozzle_expansion_ratio_stage1 double = 0

                % nozzle expansion ratio for stage2 (A_e/A_t) [-]
                nozzle_expansion_ratio_stage2 double = 0

                % density_oxidizer [kg/m^3]
                density_oxidizer double = 0

                % density_fuel [kg/m^3]
                density_fuel double = 0

                % propellant_tank_mer_oxidizer [kg/m^3]
                propellant_tank_mer_oxidizer double = 0

                % propellant_tank_mer_fuel [kg/m^3]
                propellant_tank_mer_fuel double = 0

                % propellant_tank_insulation_mer_oxidizer [kg/m^2]
                propellant_tank_insulation_mer_oxidizer double = 0

                % propellant_tank_insulation_mer_fuel [kg/m^2]
                propellant_tank_insulation_mer_fuel double = 0

                % is_solid
                is_solid logical = false
        end

        methods
                function obj = PropellantMix(name, displayname, oxidizer_mix_mass_ratio, fuel_mix_mass_ratio, specific_impulse, thrust_per_motor_stage1, thrust_per_motor_stage2, engine_exhaust_diameter_stage1, engine_exhaust_diameter_stage2, chamber_pressure_stage1, chamber_pressure_stage2, nozzle_expansion_ratio_stage1, nozzle_expansion_ratio_stage2, density_oxidizer, density_fuel, propellant_tank_mer_oxidizer, propellant_tank_mer_fuel, propellant_tank_insulation_mer_oxidizer, propellant_tank_insulation_mer_fuel, is_solid)
                        %constructor: assign all properties from inputs
                        obj.name = name;
                        obj.displayname = displayname;
                        obj.oxidizer_mix_mass_ratio = oxidizer_mix_mass_ratio;
                        obj.fuel_mix_mass_ratio = fuel_mix_mass_ratio;
                        obj.specific_impulse = specific_impulse;
                        obj.thrust_per_motor_stage1 = thrust_per_motor_stage1*1e6;
                        obj.thrust_per_motor_stage2 = thrust_per_motor_stage2*1e6;
                        obj.engine_exhaust_diameter_stage1 = engine_exhaust_diameter_stage1;
                        obj.engine_exhaust_diameter_stage2 = engine_exhaust_diameter_stage2;
                        obj.chamber_pressure_stage1 = chamber_pressure_stage1*1e6;
                        obj.chamber_pressure_stage2 = chamber_pressure_stage2*1e6;
                        obj.nozzle_expansion_ratio_stage1 = nozzle_expansion_ratio_stage1;
                        obj.nozzle_expansion_ratio_stage2 = nozzle_expansion_ratio_stage2;
                        obj.density_oxidizer = density_oxidizer;
                        obj.density_fuel = density_fuel;
                        obj.propellant_tank_mer_oxidizer = propellant_tank_mer_oxidizer;
                        obj.propellant_tank_mer_fuel = propellant_tank_mer_fuel;
                        obj.propellant_tank_insulation_mer_oxidizer = propellant_tank_insulation_mer_oxidizer;
                        obj.propellant_tank_insulation_mer_fuel = propellant_tank_insulation_mer_fuel;
                        obj.is_solid = obj.is_solid;
                end
        end
end