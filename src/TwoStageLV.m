% TwoStageLV.m
classdef TwoStageLV
        %twostagelv two-stage launch vehicle

        properties
                % propellants: propellant objects for stage1 and stage2 (1x2)
                propellants(1, 2) PropellantMix = PropellantMix("Blank", "Blank", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false)

                % xs: sweep of stage-1 delta-v fraction (1x100)
                Xs(1, 100) double = linspace(0, 1, 100)

                % ms: total stage masses (kg) for stages [1;2] at each xs (2x100)
                ms(2, 100) double = zeros(2, 100)

                % m_ins: inert mass of each stage (kg) (2x100)
                m_ins(2, 100) double = zeros(2, 100)

                % m_prs: propellant mass of each stage (kg) (2x100)
                m_prs(2, 100) double = zeros(2, 100)

                % m_0s: gross initial vehicle mass (kg) per xs (1x100)
                m_0s(1, 100) double = zeros(1, 100)

                % m_0_min: minimum gross mass found (kg)
                m_0_min double = inf

                % x_m_min: xs value that produced m_0_min
                X_m_min double = 0

                % i_m_min: index into xs for m_0_min
                i_m_min double = 1

                % costs: per-stage cost arrays (same indexing as ms) (2x100)
                costs(2, 100) double = zeros(2, 100)

                % cost_0s: total program cost per xs (1x100)
                cost_0s(1, 100) double = zeros(1, 100)

                % cost_0_min: minimum program cost found
                cost_0_min double = inf

                % x_c_min: xs that produced cost_0_min
                X_c_min double = 0

                % i_c_min: index into xs for cost_0_min
                i_c_min double = 1

                % figure handles for trends plots
                mass_fig matlab.ui.Figure
                cost_fig matlab.ui.Figure

                % struct fields are per-subsystem masses [stage1; stage2] except singletons
                m_stage_subsystem_masses struct = struct( ...
                        "payload_fairing", 0, ...
                        "propellants", [0; 0], ...
                        "propellant_tanks", [0; 0], ...
                        "propellant_tank_insulations", [0; 0], ...
                        "engines_or_casings", [0; 0], ...
                        "thrust_structures", [0; 0], ...
                        "gimbals", [0; 0], ...
                        "avionics", [0; 0], ...
                        "wirings", [0; 0], ...
                        "inter_tank_fairings", [0; 0], ...
                        "inter_stage_fairing", 0, ...
                        "aft_fairing", 0)

                m_stage_subsystem_lengths struct = struct( ...
                        "payload_fairing", 0, ...
                        "oxidizer_tanks", [0; 0], ...
                        "fuel_tanks", [0; 0], ...
                        "inter_tank_fairings", [0; 0], ...
                        "inter_stage_fairing", 0, ...
                        "aft_fairing", 0)

                % m_estimated_*: mass rollups (kg) and derived lengths (m)
                m_estimated_m_stages_margin(2, 1) double = [0; 0]   % with 30% margin on non-propellant
                m_estimated_m_total_margin double = 0
                m_estimated_m_nrec double = 0                        % payload + common fairings (no margin)
                m_estimated_m_stages(2, 1) double = [0; 0]           % no-margin per-stage
                m_estimated_m_total double = 0                        % no-margin total
                m_estimated_m_nrec_margin double = 0                 % with 30% margin on common items
                m_length_stages(2, 1) double = 0
                m_length_total double = 0
                m_necessary_engine_count double = 0
                m_estimated_c_stages(2, 1) double = [0; 0]
                m_estimated_c_total double = 0

                % c_*: (cost-oriented) mirror of the above for the second path
                c_stage_subsystem_masses struct = struct( ...
                        "payload_fairing", 0, ...
                        "propellants", [0; 0], ...
                        "propellant_tanks", [0; 0], ...
                        "propellant_tank_insulations", [0; 0], ...
                        "engines_or_casings", [0; 0], ...
                        "thrust_structures", [0; 0], ...
                        "gimbals", [0; 0], ...
                        "avionics", [0; 0], ...
                        "wirings", [0; 0], ...
                        "inter_tank_fairings", [0; 0], ...
                        "inter_stage_fairing", 0, ...
                        "aft_fairing", 0)

                c_stage_subsystem_lengths struct = struct( ...
                        "payload_fairing", 0, ...
                        "oxidizer_tanks", [0; 0], ...
                        "fuel_tanks", [0; 0], ...
                        "inter_tank_fairings", [0; 0], ...
                        "inter_stage_fairing", 0, ...
                        "aft_fairing", 0)

                c_estimated_m_stages_margin(2, 1) double = [0; 0]
                c_estimated_m_total_margin double = 0
                c_estimated_m_nrec double = 0
                c_estimated_m_stages(2, 1) double = [0; 0]
                c_estimated_m_total double = 0
                c_estimated_m_nrec_margin double = 0
                c_length_stages(2, 1) double = 0
                c_length_total double = 0
                c_necessary_engine_count double = 0
                c_estimated_c_stages(2, 1) double = [0; 0]
                c_estimated_c_total double = 0
        end

        methods
                function obj = TwoStageLV(propellants, DeltaV, m_pl, payload_reqs, delta, g, d_lv)
                        % constructor: runs part-01 trend sweep then part-02 subsystem sizing
                        arguments
                                propellants(1, 2) PropellantMix
                                DeltaV double
                                m_pl double
                                payload_reqs struct
                                delta(1, 2) double
                                g double
                                d_lv double
                        end
                        obj.propellants = propellants;

                        % part 01: delta-v split sweep, find mass/cost minima, make figs
                        obj = obj.run_part_01(DeltaV, m_pl, delta, g);

                        % part 02: engine count, subsystem masses, lengths, margins, costs
                        obj = obj.run_part_02(m_pl, payload_reqs, g, d_lv);
                end

                function obj = run_part_01(obj, DeltaV, m_pl, delta, g)
                        %run_part_01 execute part-01 mass/cost trend analysis
                        arguments
                                obj TwoStageLV
                                DeltaV double
                                m_pl double
                                delta(1, 2) double
                                g double
                        end
                        obj = obj.generate_trends(DeltaV, m_pl, delta, g);
                        obj = obj.save_mass_fig();
                        obj = obj.save_cost_fig();
                end

                function obj = run_part_02(obj, m_pl, payload_reqs, g, d_lv)
                        %run_part_02 execute part-02 sizing twice (m_* path and c_* path)
                        arguments
                                obj TwoStageLV
                                m_pl double
                                payload_reqs struct
                                g double
                                d_lv double
                        end

                        % mass-first path (m_*)
                        obj = obj.m_calculate_desired_engine_count(m_pl, payload_reqs, g, d_lv);
                        obj = obj.m_estimate_subsystem_masses(obj.m_necessary_engine_count, m_pl, payload_reqs, d_lv);
                        obj = obj.m_calculate_estimated_costs();
                        obj = obj.save_mass_table();

                        % cost-first path (c_*)
                        obj = obj.c_calculate_desired_engine_count(m_pl, payload_reqs, g, d_lv);
                        obj = obj.c_estimate_subsystem_masses(obj.c_necessary_engine_count, m_pl, payload_reqs, d_lv);
                        obj = obj.c_calculate_estimated_costs();
                        obj = obj.save_cost_table();
                end

                function obj = m_estimate_subsystem_masses(obj, n, m_pl, payload_reqs, d_lv)
                        %m_estimate_subsystem_masses compute subsystem masses for m_* path
                        % uses oxidizer/fuel split, volumes, surface areas, and curve-fits
                        arguments
                                obj TwoStageLV
                                n double
                                m_pl double
                                payload_reqs struct
                                d_lv double
                        end

                        % split prop mass into oxidizer/fuel using mix ratios
                        oxidizer_masses = [
                                obj.propellants(1).oxidizer_mix_mass_ratio / (obj.propellants(1).oxidizer_mix_mass_ratio + obj.propellants(1).fuel_mix_mass_ratio) * obj.m_prs(1, obj.i_m_min);
                                obj.propellants(2).oxidizer_mix_mass_ratio / (obj.propellants(2).oxidizer_mix_mass_ratio + obj.propellants(2).fuel_mix_mass_ratio) * obj.m_prs(2, obj.i_m_min);
                        ];

                        fuel_masses = [
                                obj.propellants(1).fuel_mix_mass_ratio / (obj.propellants(1).oxidizer_mix_mass_ratio + obj.propellants(1).fuel_mix_mass_ratio) * obj.m_prs(1, obj.i_m_min);
                                obj.propellants(2).fuel_mix_mass_ratio / (obj.propellants(2).oxidizer_mix_mass_ratio + obj.propellants(2).fuel_mix_mass_ratio) * obj.m_prs(2, obj.i_m_min);
                        ];

                        % convert masses to volumes via densities
                        oxidizer_volumes = [
                                oxidizer_masses(1)/obj.propellants(1).density_oxidizer;
                                oxidizer_masses(2)/obj.propellants(2).density_oxidizer;
                        ];
                        fuel_volumes = [
                                fuel_masses(1)/obj.propellants(1).density_fuel;
                                fuel_masses(2)/obj.propellants(2).density_fuel;
                        ];

                        % geometric assumptions for tanks and fairings
                        r_tank = d_lv/2; % desired tank radius is half of lv diameter
                        r_payload = payload_reqs.diameter/2;

                        % payload cylinder + cone + hemisphere (area/height model)
                        frustum_height = 5/2*(r_tank - r_payload);
                        payload_height = frustum_height + payload_reqs.height + r_payload;
                        payload_area = payload_reqs.height*2*pi*r_tank + (pi*r_tank*sqrt((5/2*r_tank)^2 + r_tank^2) - pi*r_payload*sqrt((5/2*r_payload)^2 + r_payload^2)) + 2*pi*r_payload^2;

                        % cylindrical tank heights (vol = pi*r^2*h) and lateral areas (sa = 2*pi*r*h)
                        oxidizer_tank_heights = oxidizer_volumes./(pi*r_tank^2);
                        oxidizer_tank_areas = 2*pi*r_tank*oxidizer_tank_heights;

                        fuel_tank_heights = fuel_volumes./(pi*r_tank^2);
                        fuel_tank_areas = 2*pi*r_tank*fuel_tank_heights;

                        % spacers and interstage regions
                        inter_tank_heights = [0.5*r_tank; 0.5*r_tank];
                        inter_tank_areas = 2*pi*r_tank*inter_tank_heights;

                        inter_stage_height = 0.5*r_tank;
                        inter_stage_area = 2*pi*r_tank*inter_stage_height;

                        aft_height = 0.5*r_tank;
                        aft_area = 2*pi*r_tank*aft_height;

                        % length bookkeeping (per-stage and total)
                        obj.m_stage_subsystem_lengths.payload_fairing = payload_height;
                        obj.m_stage_subsystem_lengths.oxidizer_tanks = oxidizer_tank_heights;
                        obj.m_stage_subsystem_lengths.fuel_tanks = fuel_tank_heights;
                        obj.m_stage_subsystem_lengths.inter_tank_fairings = inter_tank_heights;
                        obj.m_stage_subsystem_lengths.inter_stage_fairing = inter_stage_height;
                        obj.m_stage_subsystem_lengths.aft_fairing = aft_height;

                        obj.m_length_stages = oxidizer_tank_heights+fuel_tank_heights+inter_tank_heights;
                        obj.m_length_total = payload_height+sum(obj.m_length_stages)+inter_stage_height+aft_height;

                        % subsystem masses via fits/heuristics
                        obj.m_stage_subsystem_masses.propellants = [obj.m_prs(1, obj.i_m_min); obj.m_prs(2, obj.i_m_min)];
                        obj.m_stage_subsystem_masses.propellant_tanks = [
                                obj.propellants(1).propellant_tank_mer_oxidizer*oxidizer_volumes(1) + obj.propellants(1).propellant_tank_mer_fuel*fuel_volumes(1);
                                obj.propellants(2).propellant_tank_mer_oxidizer*oxidizer_volumes(2) + obj.propellants(2).propellant_tank_mer_fuel*fuel_volumes(2);
                        ];
                        obj.m_stage_subsystem_masses.propellant_tank_insulations = [
                                obj.propellants(1).propellant_tank_insulation_mer_oxidizer*oxidizer_tank_areas(1) + obj.propellants(1).propellant_tank_insulation_mer_fuel*fuel_tank_areas(1);
                                obj.propellants(2).propellant_tank_insulation_mer_oxidizer*oxidizer_tank_areas(2) + obj.propellants(2).propellant_tank_insulation_mer_fuel*fuel_tank_areas(2);
                        ];
                        obj.m_stage_subsystem_masses.engines_or_casings = [
                                n*m_calculate_stage_engine_or_casing_mass(obj, 1, obj.propellants(1).thrust_per_motor_stage1, obj.propellants(1).nozzle_expansion_ratio_stage1);
                                n*m_calculate_stage_engine_or_casing_mass(obj, 2, obj.propellants(2).thrust_per_motor_stage2, obj.propellants(2).nozzle_expansion_ratio_stage2);
                        ];
                        obj.m_stage_subsystem_masses.thrust_structures = [
                                n*2.55e-4*obj.propellants(1).thrust_per_motor_stage1;
                                n*2.55e-4*obj.propellants(2).thrust_per_motor_stage2;
                        ];
                        obj.m_stage_subsystem_masses.gimbals = [
                                n*237.8*(obj.propellants(1).thrust_per_motor_stage1/obj.propellants(1).chamber_pressure_stage1)^0.9375;
                                n*237.8*(obj.propellants(2).thrust_per_motor_stage2/obj.propellants(2).chamber_pressure_stage2)^0.9375;
                        ];
                        obj.m_stage_subsystem_masses.avionics = [
                                10*obj.ms(1, obj.i_m_min)^0.361;
                                10*obj.ms(2, obj.i_m_min)^0.361;
                        ];
                        obj.m_stage_subsystem_masses.wirings = [
                                1.058*sqrt(obj.ms(1, obj.i_m_min))*obj.m_length_stages(1)^(0.25);
                                1.058*sqrt(obj.ms(2, obj.i_m_min))*obj.m_length_stages(2)^(0.25);
                        ];
                        obj.m_stage_subsystem_masses.payload_fairing = 4.95*payload_area^1.15;
                        obj.m_stage_subsystem_masses.inter_tank_fairings = [
                                4.95*inter_tank_areas(1)^1.15;
                                4.95*inter_tank_areas(2)^1.15;
                        ];
                        obj.m_stage_subsystem_masses.inter_stage_fairing = 4.95*inter_stage_area^1.15;
                        obj.m_stage_subsystem_masses.aft_fairing = 4.95*aft_area^1.15;

                        % rollups (no-margin and with 30% margin applied to non-propellant)
                        obj.m_estimated_m_stages = [
                                obj.m_stage_subsystem_masses.propellants(1)+obj.m_stage_subsystem_masses.propellant_tanks(1)+obj.m_stage_subsystem_masses.propellant_tank_insulations(1)+obj.m_stage_subsystem_masses.engines_or_casings(1)+obj.m_stage_subsystem_masses.thrust_structures(1)+obj.m_stage_subsystem_masses.gimbals(1)+obj.m_stage_subsystem_masses.avionics(1)+obj.m_stage_subsystem_masses.wirings(1)+obj.m_stage_subsystem_masses.inter_tank_fairings(1);
                                obj.m_stage_subsystem_masses.propellants(2)+obj.m_stage_subsystem_masses.propellant_tanks(2)+obj.m_stage_subsystem_masses.propellant_tank_insulations(2)+obj.m_stage_subsystem_masses.engines_or_casings(2)+obj.m_stage_subsystem_masses.thrust_structures(2)+obj.m_stage_subsystem_masses.gimbals(2)+obj.m_stage_subsystem_masses.avionics(2)+obj.m_stage_subsystem_masses.wirings(2)+obj.m_stage_subsystem_masses.inter_tank_fairings(2);
                        ];
                        obj.m_estimated_m_nrec = m_pl + obj.m_stage_subsystem_masses.payload_fairing+obj.m_stage_subsystem_masses.inter_stage_fairing+obj.m_stage_subsystem_masses.aft_fairing;
                        obj.m_estimated_m_total = sum(obj.m_estimated_m_stages)+obj.m_estimated_m_nrec;

                        obj.m_estimated_m_stages_margin = [
                                obj.m_stage_subsystem_masses.propellants(1) + 1.3*(obj.m_stage_subsystem_masses.propellant_tanks(1)+obj.m_stage_subsystem_masses.propellant_tank_insulations(1)+obj.m_stage_subsystem_masses.engines_or_casings(1)+obj.m_stage_subsystem_masses.thrust_structures(1)+obj.m_stage_subsystem_masses.gimbals(1)+obj.m_stage_subsystem_masses.avionics(1)+obj.m_stage_subsystem_masses.wirings(1)+obj.m_stage_subsystem_masses.inter_tank_fairings(1));
                                obj.m_stage_subsystem_masses.propellants(2) + 1.3*(obj.m_stage_subsystem_masses.propellant_tanks(2)+obj.m_stage_subsystem_masses.propellant_tank_insulations(2)+obj.m_stage_subsystem_masses.engines_or_casings(2)+obj.m_stage_subsystem_masses.thrust_structures(2)+obj.m_stage_subsystem_masses.gimbals(2)+obj.m_stage_subsystem_masses.avionics(2)+obj.m_stage_subsystem_masses.wirings(2)+obj.m_stage_subsystem_masses.inter_tank_fairings(2));
                        ];
                        obj.m_estimated_m_nrec_margin = m_pl + 1.3*(obj.m_stage_subsystem_masses.payload_fairing+obj.m_stage_subsystem_masses.inter_stage_fairing+obj.m_stage_subsystem_masses.aft_fairing);
                        obj.m_estimated_m_total_margin = sum(obj.m_estimated_m_stages_margin)+obj.m_estimated_m_nrec_margin;
                end

                function m_e_or_c = m_calculate_stage_engine_or_casing_mass(obj, stage, T_N, NER)
                        %m_calculate_stage_engine_or_casing_mass liquid engine or solid case mass
                        % returns casing mass ~ 13.5% of propellant for solids; otherwise engine fit
                        arguments
                                obj TwoStageLV
                                stage double
                                T_N double
                                NER double
                        end
                        if obj.propellants(stage).is_solid
                                m_e_or_c = 0.135*obj.m_stage_subsystem_masses.propellants(stage);
                        else
                                m_e_or_c = 7.81e-4*T_N + 3.37e-5*T_N*sqrt(NER) + 59;
                        end
                end

                function obj = m_calculate_desired_engine_count(obj, m_pl, payload_reqs, g, d_lv)
                        %m_calculate_desired_engine_count find smallest n meeting thrust margins
                        % stage 1: t/w >= 1.3 on 1.3*total_with_margin; stage 2: t/w >= 0.76 on upper mass
                        arguments
                                obj TwoStageLV
                                m_pl double
                                payload_reqs struct
                                g double
                                d_lv double
                        end

                        n = 1;
                        valid_engine_count = false;
                        while ~valid_engine_count
                                obj = obj.m_estimate_subsystem_masses(n, m_pl, payload_reqs, d_lv);

                                T_valid = [
                                        obj.propellants(1).thrust_per_motor_stage1 * n >= 1.3 * (1.3*obj.m_estimated_m_total_margin) * g;
                                        obj.propellants(2).thrust_per_motor_stage2 * n >= 0.76 * (1.3*(obj.m_estimated_m_stages_margin(2)+obj.m_estimated_m_nrec_margin)) * g;
                                ];

                                valid_engine_count = all(T_valid);
                                if ~valid_engine_count
                                        n = n + 1;
                                end
                        end

                        obj.m_necessary_engine_count = n;
                end

                function obj = m_calculate_estimated_costs(obj)
                        %m_calculate_estimated_costs cost via power law on mass-with-margin
                        arguments
                                obj TwoStageLV
                        end
                        
                        obj.m_estimated_c_stages = [
                                13.52*(obj.m_estimated_m_stages_margin(1)-obj.m_stage_subsystem_masses.propellants(1))^(0.55)*1e6;
                                13.52*(obj.m_estimated_m_stages_margin(2)-obj.m_stage_subsystem_masses.propellants(2))^(0.55)*1e6;
                        ];
                        obj.m_estimated_c_total = 13.52*(obj.m_estimated_m_total_margin-sum(obj.m_stage_subsystem_masses.propellants))^(0.55)*1e6;
                end

                function obj = save_mass_table(obj)
                        %save_mass_table write csv for the m_* path breakdown
                        mass_table = table( ...
                                ["-"; "-"; obj.m_stage_subsystem_masses.payload_fairing], ...
                                ["-"; "-"; obj.m_stage_subsystem_lengths.payload_fairing], ...
                                [obj.m_stage_subsystem_masses.propellants(1); obj.m_stage_subsystem_masses.propellants(2); sum(obj.m_stage_subsystem_masses.propellants)], ...
                                [obj.m_stage_subsystem_masses.propellant_tanks(1); obj.m_stage_subsystem_masses.propellant_tanks(2); sum(obj.m_stage_subsystem_masses.propellant_tanks)], ...
                                [obj.m_stage_subsystem_masses.propellant_tank_insulations(1); obj.m_stage_subsystem_masses.propellant_tank_insulations(2); sum(obj.m_stage_subsystem_masses.propellant_tank_insulations)], ...
                                [obj.m_stage_subsystem_lengths.oxidizer_tanks(1); obj.m_stage_subsystem_lengths.oxidizer_tanks(2); sum(obj.m_stage_subsystem_lengths.oxidizer_tanks)], ...
                                [obj.m_stage_subsystem_lengths.fuel_tanks(1); obj.m_stage_subsystem_lengths.fuel_tanks(2); sum(obj.m_stage_subsystem_lengths.fuel_tanks)], ...
                                [obj.m_stage_subsystem_masses.engines_or_casings(1); obj.m_stage_subsystem_masses.engines_or_casings(2); sum(obj.m_stage_subsystem_masses.engines_or_casings)], ...
                                [obj.m_stage_subsystem_masses.thrust_structures(1); obj.m_stage_subsystem_masses.thrust_structures(2); sum(obj.m_stage_subsystem_masses.thrust_structures)], ...
                                [obj.m_stage_subsystem_masses.gimbals(1); obj.m_stage_subsystem_masses.gimbals(2); sum(obj.m_stage_subsystem_masses.gimbals)], ...
                                [obj.m_stage_subsystem_masses.avionics(1); obj.m_stage_subsystem_masses.avionics(2); sum(obj.m_stage_subsystem_masses.avionics)], ...
                                [obj.m_stage_subsystem_masses.wirings(1); obj.m_stage_subsystem_masses.wirings(2); sum(obj.m_stage_subsystem_masses.wirings)], ...
                                [obj.m_stage_subsystem_masses.inter_tank_fairings(1); obj.m_stage_subsystem_masses.inter_tank_fairings(2); sum(obj.m_stage_subsystem_masses.inter_tank_fairings)], ...
                                [obj.m_stage_subsystem_lengths.inter_tank_fairings(1); obj.m_stage_subsystem_lengths.inter_tank_fairings(2); sum(obj.m_stage_subsystem_lengths.inter_tank_fairings)], ...
                                ["-"; "-"; obj.m_stage_subsystem_masses.inter_stage_fairing], ...
                                ["-"; "-"; obj.m_stage_subsystem_lengths.inter_stage_fairing], ...
                                ["-"; "-"; obj.m_stage_subsystem_masses.aft_fairing], ...
                                ["-"; "-"; obj.m_stage_subsystem_lengths.aft_fairing], ...
                                [obj.m_estimated_m_stages(1); obj.m_estimated_m_stages(2); obj.m_estimated_m_total], ...
                                [obj.m_estimated_m_stages_margin(1); obj.m_estimated_m_stages_margin(2); obj.m_estimated_m_total_margin], ...
                                [obj.m_estimated_m_stages(1)-obj.m_stage_subsystem_masses.propellants(1); obj.m_estimated_m_stages(2)-obj.m_stage_subsystem_masses.propellants(2); obj.m_estimated_m_total-sum(obj.m_stage_subsystem_masses.propellants)], ...
                                [obj.m_estimated_m_stages_margin(1)-obj.m_stage_subsystem_masses.propellants(1); obj.m_estimated_m_stages_margin(2)-obj.m_stage_subsystem_masses.propellants(2); obj.m_estimated_m_total_margin-sum(obj.m_stage_subsystem_masses.propellants)], ...
                                [obj.m_length_stages(1); obj.m_length_stages(2); obj.m_length_total], ...
                                ["-"; "-"; obj.m_necessary_engine_count], ...
                                [obj.m_estimated_c_stages(1)/1e9; obj.m_estimated_c_stages(2)/1e9; obj.m_estimated_c_total/1e9], ...
                                RowNames=["Stage 1", "Stage 2", "Total"], ...
                                VariableNames=[ ...
                                        "Payload Fairing Mass [kg]";
                                        "Payload Fairing Length [m]";
                                        "Propellant Mass [kg]";
                                        "Propellant Tank Mass [kg]";
                                        "Propellant Tank Insulation Mass [kg]";
                                        "Oxidizer Tank Length [m]";
                                        "Fuel Tank Length [m]";
                                        "Engines/Casing Mass [kg]";
                                        "Thrust Structure Mass [kg]";
                                        "Gimbal Mass [kg]";
                                        "Avionic Mass [kg]";
                                        "Wiring Mass [kg]";
                                        "Inter-Tank Fairing Mass [kg]";
                                        "Inter-Tank Fairing Length [m]";
                                        "Inter-Stage Fairing Mass [kg]";
                                        "Inter-Stage Fairing Length [m]";
                                        "Aft Fairing Mass [kg]";
                                        "Aft Fairing Length [m]";
                                        "Stage Mass [kg]";
                                        "Stage Mass w/ Margin [kg]";
                                        "Inert Mass [kg]";
                                        "Inert Mass w/ Margin [kg]";
                                        "Stage Length [m]";
                                        "Engine Count [-]";
                                        "Stage Cost [B$USD]";
                                ]);

                        % format numeric columns to 3 significant digits
                        idx = vartype("numeric");
                        mass_table{:, idx} = round(mass_table{:, idx}, 3, "significant");

                        % write csv under ./output/<s2>/mass/s1 <s1> - s2 <s2>.csv
                        filename = sprintf("./output/%2$s/mass/s1 %1$s - s2 %2$s.csv", obj.propellants(1).name, obj.propellants(2).name);
                        writetable(mass_table, filename, WriteRowNames=true, WriteVariableNames=true);
                end

                function obj = c_estimate_subsystem_masses(obj, n, m_pl, payload_reqs, d_lv)
                        %c_estimate_subsystem_masses compute subsystem masses for c_* path
                        arguments
                                obj TwoStageLV
                                n double
                                m_pl double
                                payload_reqs struct
                                d_lv double
                        end

                        % same flow as m_* path, using separate c_* fields
                        oxidizer_masses = [
                                obj.propellants(1).oxidizer_mix_mass_ratio / (obj.propellants(1).oxidizer_mix_mass_ratio + obj.propellants(1).fuel_mix_mass_ratio) * obj.m_prs(1, obj.i_c_min);
                                obj.propellants(2).oxidizer_mix_mass_ratio / (obj.propellants(2).oxidizer_mix_mass_ratio + obj.propellants(2).fuel_mix_mass_ratio) * obj.m_prs(2, obj.i_c_min);
                        ];
                        fuel_masses = [
                                obj.propellants(1).fuel_mix_mass_ratio / (obj.propellants(1).oxidizer_mix_mass_ratio + obj.propellants(1).fuel_mix_mass_ratio) * obj.m_prs(1, obj.i_c_min);
                                obj.propellants(2).fuel_mix_mass_ratio / (obj.propellants(2).oxidizer_mix_mass_ratio + obj.propellants(2).fuel_mix_mass_ratio) * obj.m_prs(2, obj.i_c_min);
                        ];
                        oxidizer_volumes = [
                                oxidizer_masses(1)/obj.propellants(1).density_oxidizer;
                                oxidizer_masses(2)/obj.propellants(2).density_oxidizer;
                        ];
                        fuel_volumes = [
                                fuel_masses(1)/obj.propellants(1).density_fuel;
                                fuel_masses(2)/obj.propellants(2).density_fuel;
                        ];


                        % same geometry assumptions
                        r_tank = d_lv/2;
                        r_payload = payload_reqs.diameter/2;

                        frustum_height = 5/2*(r_tank - r_payload);
                        payload_height = frustum_height + payload_reqs.height + r_payload;
                        payload_area = payload_reqs.height*2*pi*r_tank + (pi*r_tank*sqrt((5/2*r_tank)^2 + r_tank^2) - pi*r_payload*sqrt((5/2*r_payload)^2 + r_payload^2)) + 2*pi*r_payload^2;

                        oxidizer_tank_heights = oxidizer_volumes./(pi*r_tank^2);
                        oxidizer_tank_areas = 2*pi*r_tank*oxidizer_tank_heights;

                        fuel_tank_heights = fuel_volumes./(pi*r_tank^2);
                        fuel_tank_areas = 2*pi*r_tank*fuel_tank_heights;

                        inter_tank_heights = [0.5*r_tank; 0.5*r_tank];
                        inter_tank_areas = 2*pi*r_tank*inter_tank_heights;

                        inter_stage_height = 0.5*r_tank;
                        inter_stage_area = 2*pi*r_tank*inter_stage_height;

                        aft_height = 0.5*r_tank;
                        aft_area = 2*pi*r_tank*aft_height;

                        obj.c_stage_subsystem_lengths.payload_fairing = payload_height;
                        obj.c_stage_subsystem_lengths.oxidizer_tanks = oxidizer_tank_heights;
                        obj.c_stage_subsystem_lengths.fuel_tanks = fuel_tank_heights;
                        obj.c_stage_subsystem_lengths.inter_tank_fairings = inter_tank_heights;
                        obj.c_stage_subsystem_lengths.inter_stage_fairing = inter_stage_height;
                        obj.c_stage_subsystem_lengths.aft_fairing = aft_height;

                        obj.c_length_stages = oxidizer_tank_heights+fuel_tank_heights+inter_tank_heights;
                        obj.c_length_total = payload_height+sum(obj.c_length_stages)+inter_stage_height+aft_height;

                        obj.c_stage_subsystem_masses.propellants = [obj.m_prs(1, obj.i_c_min); obj.m_prs(2, obj.i_c_min)];
                        obj.c_stage_subsystem_masses.propellant_tanks = [
                                obj.propellants(1).propellant_tank_mer_oxidizer*oxidizer_volumes(1) + obj.propellants(1).propellant_tank_mer_fuel*fuel_volumes(1);
                                obj.propellants(2).propellant_tank_mer_oxidizer*oxidizer_volumes(2) + obj.propellants(2).propellant_tank_mer_fuel*fuel_volumes(2);
                        ];
                        obj.c_stage_subsystem_masses.propellant_tank_insulations = [
                                obj.propellants(1).propellant_tank_insulation_mer_oxidizer*oxidizer_tank_areas(1) + obj.propellants(1).propellant_tank_insulation_mer_fuel*fuel_tank_areas(1);
                                obj.propellants(2).propellant_tank_insulation_mer_oxidizer*oxidizer_tank_areas(2) + obj.propellants(2).propellant_tank_insulation_mer_fuel*fuel_tank_areas(2);
                        ];
                        obj.c_stage_subsystem_masses.engines_or_casings = [
                                n*c_calculate_stage_engine_or_casing_mass(obj, 1, obj.propellants(1).thrust_per_motor_stage1, obj.propellants(1).nozzle_expansion_ratio_stage1);
                                n*c_calculate_stage_engine_or_casing_mass(obj, 2, obj.propellants(2).thrust_per_motor_stage2, obj.propellants(2).nozzle_expansion_ratio_stage2);
                        ];
                        obj.c_stage_subsystem_masses.thrust_structures = [
                                n*2.55e-4*obj.propellants(1).thrust_per_motor_stage1;
                                n*2.55e-4*obj.propellants(2).thrust_per_motor_stage2;
                        ];
                        obj.c_stage_subsystem_masses.gimbals = [
                                n*237.8*(obj.propellants(1).thrust_per_motor_stage1/obj.propellants(1).chamber_pressure_stage1)^0.9375;
                                n*237.8*(obj.propellants(2).thrust_per_motor_stage2/obj.propellants(2).chamber_pressure_stage2)^0.9375;
                        ];
                        obj.c_stage_subsystem_masses.avionics = [
                                10*obj.ms(1, obj.i_c_min)^0.361;
                                10*obj.ms(2, obj.i_c_min)^0.361;
                        ];
                        obj.c_stage_subsystem_masses.wirings = [
                                1.058*sqrt(obj.ms(1, obj.i_c_min))*obj.c_length_stages(1)^(0.25);
                                1.058*sqrt(obj.ms(2, obj.i_c_min))*obj.c_length_stages(2)^(0.25);
                        ];
                        obj.c_stage_subsystem_masses.payload_fairing = 4.95*payload_area^1.15;
                        obj.c_stage_subsystem_masses.inter_tank_fairings = [
                                4.95*inter_tank_areas(1)^1.15;
                                4.95*inter_tank_areas(2)^1.15;
                        ];
                        obj.c_stage_subsystem_masses.inter_stage_fairing = 4.95*inter_stage_area^1.15;
                        obj.c_stage_subsystem_masses.aft_fairing = 4.95*aft_area^1.15;

                        obj.c_estimated_m_stages = [
                                obj.c_stage_subsystem_masses.propellants(1)+obj.c_stage_subsystem_masses.propellant_tanks(1)+obj.c_stage_subsystem_masses.propellant_tank_insulations(1)+obj.c_stage_subsystem_masses.engines_or_casings(1)+obj.c_stage_subsystem_masses.thrust_structures(1)+obj.c_stage_subsystem_masses.gimbals(1)+obj.c_stage_subsystem_masses.avionics(1)+obj.c_stage_subsystem_masses.wirings(1)+obj.c_stage_subsystem_masses.inter_tank_fairings(1);
                                obj.c_stage_subsystem_masses.propellants(2)+obj.c_stage_subsystem_masses.propellant_tanks(2)+obj.c_stage_subsystem_masses.propellant_tank_insulations(2)+obj.c_stage_subsystem_masses.engines_or_casings(2)+obj.c_stage_subsystem_masses.thrust_structures(2)+obj.c_stage_subsystem_masses.gimbals(2)+obj.c_stage_subsystem_masses.avionics(2)+obj.c_stage_subsystem_masses.wirings(2)+obj.c_stage_subsystem_masses.inter_tank_fairings(2);
                        ];
                        obj.c_estimated_m_nrec = m_pl + obj.c_stage_subsystem_masses.payload_fairing+obj.c_stage_subsystem_masses.inter_stage_fairing+obj.c_stage_subsystem_masses.aft_fairing;
                        obj.c_estimated_m_total = sum(obj.c_estimated_m_stages)+obj.c_estimated_m_nrec;

                        obj.c_estimated_m_stages_margin = [
                                obj.c_stage_subsystem_masses.propellants(1) + 1.3*(obj.c_stage_subsystem_masses.propellant_tanks(1)+obj.c_stage_subsystem_masses.propellant_tank_insulations(1)+obj.c_stage_subsystem_masses.engines_or_casings(1)+obj.c_stage_subsystem_masses.thrust_structures(1)+obj.c_stage_subsystem_masses.gimbals(1)+obj.c_stage_subsystem_masses.avionics(1)+obj.c_stage_subsystem_masses.wirings(1)+obj.c_stage_subsystem_masses.inter_tank_fairings(1));
                                obj.c_stage_subsystem_masses.propellants(2) + 1.3*(obj.c_stage_subsystem_masses.propellant_tanks(2)+obj.c_stage_subsystem_masses.propellant_tank_insulations(2)+obj.c_stage_subsystem_masses.engines_or_casings(2)+obj.c_stage_subsystem_masses.thrust_structures(2)+obj.c_stage_subsystem_masses.gimbals(2)+obj.c_stage_subsystem_masses.avionics(2)+obj.c_stage_subsystem_masses.wirings(2)+obj.c_stage_subsystem_masses.inter_tank_fairings(2));
                        ];
                        obj.c_estimated_m_nrec_margin = m_pl + 1.3*(obj.c_stage_subsystem_masses.payload_fairing+obj.c_stage_subsystem_masses.inter_stage_fairing+obj.c_stage_subsystem_masses.aft_fairing);
                        obj.c_estimated_m_total_margin = sum(obj.c_estimated_m_stages_margin)+obj.c_estimated_m_nrec_margin;
                end

                function m_e_or_c = c_calculate_stage_engine_or_casing_mass(obj, stage, T_N, NER)
                        %c_calculate_stage_engine_or_casing_mass identical to m_* but for c_* path
                        arguments
                                obj TwoStageLV
                                stage double
                                T_N double
                                NER double
                        end
                        if obj.propellants(stage).is_solid
                                m_e_or_c = 0.135*obj.c_stage_subsystem_masses.propellants(stage);
                        else
                                m_e_or_c = 7.81e-4*T_N + 3.37e-5*T_N*sqrt(NER) + 59;
                        end
                end

                function obj = c_calculate_desired_engine_count(obj, m_pl, payload_reqs, g, d_lv)
                        %c_calculate_desired_engine_count same thrust-margin search for c_* path
                        arguments
                                obj TwoStageLV
                                m_pl double
                                payload_reqs struct
                                g double
                                d_lv double
                        end

                        n = 1;
                        valid_engine_count = false;
                        while ~valid_engine_count
                                obj = obj.c_estimate_subsystem_masses(n, m_pl, payload_reqs, d_lv);

                                T_valid = [
                                        obj.propellants(1).thrust_per_motor_stage1 * n >= 1.3 * (1.3*obj.c_estimated_m_total_margin) * g;
                                        obj.propellants(2).thrust_per_motor_stage2 * n >= 0.76 * (1.3*(obj.c_estimated_m_stages_margin(2)+obj.c_estimated_m_nrec_margin)) * g;
                                ];

                                valid_engine_count = all(T_valid);
                                if ~valid_engine_count
                                        n = n + 1;
                                end
                        end

                        obj.c_necessary_engine_count = n;
                end

                function obj = c_calculate_estimated_costs(obj)
                        %c_calculate_estimated_costs cost for c_* path based on masses-with-margin
                        arguments
                                obj TwoStageLV
                        end
                        obj.c_estimated_c_stages = [
                                13.52*(obj.c_estimated_m_stages_margin(1)-obj.c_stage_subsystem_masses.propellants(1))^(0.55)*1e6;
                                13.52*(obj.c_estimated_m_stages_margin(2)-obj.c_stage_subsystem_masses.propellants(2))^(0.55)*1e6;
                        ];
                        obj.c_estimated_c_total = 13.52*(obj.c_estimated_m_total_margin-sum(obj.c_stage_subsystem_masses.propellants))^(0.55)*1e6;
                end

                function obj = save_cost_table(obj)
                        %save_cost_table write csv for the c_* path breakdown
                        cost_table = table( ...
                                ["-"; "-"; obj.c_stage_subsystem_masses.payload_fairing], ...
                                ["-"; "-"; obj.c_stage_subsystem_lengths.payload_fairing], ...
                                [obj.c_stage_subsystem_masses.propellants(1); obj.c_stage_subsystem_masses.propellants(2); sum(obj.c_stage_subsystem_masses.propellants)], ...
                                [obj.c_stage_subsystem_masses.propellant_tanks(1); obj.c_stage_subsystem_masses.propellant_tanks(2); sum(obj.c_stage_subsystem_masses.propellant_tanks)], ...
                                [obj.c_stage_subsystem_masses.propellant_tank_insulations(1); obj.c_stage_subsystem_masses.propellant_tank_insulations(2); sum(obj.c_stage_subsystem_masses.propellant_tank_insulations)], ...
                                [obj.c_stage_subsystem_lengths.oxidizer_tanks(1); obj.c_stage_subsystem_lengths.oxidizer_tanks(2); sum(obj.c_stage_subsystem_lengths.oxidizer_tanks)], ...
                                [obj.c_stage_subsystem_lengths.fuel_tanks(1); obj.c_stage_subsystem_lengths.fuel_tanks(2); sum(obj.c_stage_subsystem_lengths.fuel_tanks)], ...
                                [obj.c_stage_subsystem_masses.engines_or_casings(1); obj.c_stage_subsystem_masses.engines_or_casings(2); sum(obj.c_stage_subsystem_masses.engines_or_casings)], ...
                                [obj.c_stage_subsystem_masses.thrust_structures(1); obj.c_stage_subsystem_masses.thrust_structures(2); sum(obj.c_stage_subsystem_masses.thrust_structures)], ...
                                [obj.c_stage_subsystem_masses.gimbals(1); obj.c_stage_subsystem_masses.gimbals(2); sum(obj.c_stage_subsystem_masses.gimbals)], ...
                                [obj.c_stage_subsystem_masses.avionics(1); obj.c_stage_subsystem_masses.avionics(2); sum(obj.c_stage_subsystem_masses.avionics)], ...
                                [obj.c_stage_subsystem_masses.wirings(1); obj.c_stage_subsystem_masses.wirings(2); sum(obj.c_stage_subsystem_masses.wirings)], ...
                                [obj.c_stage_subsystem_masses.inter_tank_fairings(1); obj.c_stage_subsystem_masses.inter_tank_fairings(2); sum(obj.c_stage_subsystem_masses.inter_tank_fairings)], ...
                                [obj.c_stage_subsystem_lengths.inter_tank_fairings(1); obj.c_stage_subsystem_lengths.inter_tank_fairings(2); sum(obj.c_stage_subsystem_lengths.inter_tank_fairings)], ...
                                ["-"; "-"; obj.c_stage_subsystem_masses.inter_stage_fairing], ...
                                ["-"; "-"; obj.c_stage_subsystem_lengths.inter_stage_fairing], ...
                                ["-"; "-"; obj.c_stage_subsystem_masses.aft_fairing], ...
                                ["-"; "-"; obj.c_stage_subsystem_lengths.aft_fairing], ...
                                [obj.c_estimated_m_stages(1); obj.c_estimated_m_stages(2); obj.c_estimated_m_total], ...
                                [obj.c_estimated_m_stages_margin(1); obj.c_estimated_m_stages_margin(2); obj.c_estimated_m_total_margin], ...
                                [obj.c_estimated_m_stages(1)-obj.c_stage_subsystem_masses.propellants(1); obj.c_estimated_m_stages(2)-obj.c_stage_subsystem_masses.propellants(2); obj.c_estimated_m_total-sum(obj.c_stage_subsystem_masses.propellants)], ...
                                [obj.c_estimated_m_stages_margin(1)-obj.c_stage_subsystem_masses.propellants(1); obj.c_estimated_m_stages_margin(2)-obj.c_stage_subsystem_masses.propellants(2); obj.c_estimated_m_total_margin-sum(obj.c_stage_subsystem_masses.propellants)], ...
                                [obj.c_length_stages(1); obj.c_length_stages(2); obj.c_length_total], ...
                                ["-"; "-"; obj.c_necessary_engine_count], ...
                                [obj.c_estimated_c_stages(1)/1e9; obj.c_estimated_c_stages(2)/1e9; obj.c_estimated_c_total/1e9], ...
                                RowNames=["Stage 1", "Stage 2", "Total"], ...
                                VariableNames=[ ...
                                        "Payload Fairing Mass [kg]";
                                        "Payload Fairing Length [m]";
                                        "Propellant Mass [kg]";
                                        "Propellant Tank Mass [kg]";
                                        "Propellant Tank Insulation Mass [kg]";
                                        "Oxidizer Tank Length [m]";
                                        "Fuel Tank Length [m]";
                                        "Engines/Casing Mass [kg]";
                                        "Thrust Structure Mass [kg]";
                                        "Gimbal Mass [kg]";
                                        "Avionic Mass [kg]";
                                        "Wiring Mass [kg]";
                                        "Inter-Tank Fairing Mass [kg]";
                                        "Inter-Tank Fairing Length [m]";
                                        "Inter-Stage Fairing Mass [kg]";
                                        "Inter-Stage Fairing Length [m]";
                                        "Aft Fairing Mass [kg]";
                                        "Aft Fairing Length [m]";
                                        "Stage Mass [kg]";
                                        "Stage Mass w/ Margin [kg]";
                                        "Inert Mass [kg]";
                                        "Inert Mass w/ Margin [kg]";
                                        "Stage Length [m]";
                                        "Engine Count [-]";
                                        "Stage Cost [B$USD]";
                                ]);

                        % format numeric columns to 3 significant digits
                        idx = vartype("numeric");
                        cost_table{:, idx} = round(cost_table{:, idx}, 3, "significant");

                        % write csv under ./output/<s2>/cost/s1 <s1> - s2 <s2>.csv
                        filename = sprintf("./output/%2$s/cost/s1 %1$s - s2 %2$s.csv", obj.propellants(1).name, obj.propellants(2).name);
                        writetable(cost_table, filename, WriteRowNames=true, WriteVariableNames=true);
                end

                function [m, m_in, m_pr, m_0] = calculate_stage_masses(obj, X, DeltaV, m_pl, delta, g)
                        %calculate_stage_masses compute per-stage masses for a given x split
                        arguments
                                obj TwoStageLV
                                X double
                                DeltaV double
                                m_pl double
                                delta(1, 2) double
                                g double
                        end

                        % ve per stage (m/s): ve = isp * g
                        Ve(:) = [g*obj.propellants(1).specific_impulse, g*obj.propellants(2).specific_impulse];

                        % delta-v split across stages
                        DV(:) = [X*DeltaV, (1-X)*DeltaV];

                        % rocket eqn mass ratios
                        r(:) = [exp(-DV(1)/Ve(1)), exp(-DV(2)/Ve(2))];

                        % helper lambda = r - delta (closed-form solution)
                        lambda(:) = [r(1) - delta(1), r(2) - delta(2)];

                        % initial masses per stage (propagate upward from payload)
                        m_o(2) = m_pl/lambda(2);
                        m_o(1) = m_o(2)/lambda(1);

                        % inert per stage
                        m_in(:) = [delta(1)*m_o(1), delta(2)*m_o(2)];

                        % propellant per stage
                        m_pr(:) = [m_o(1) - m_o(2) - m_in(1), m_o(2) - m_pl - m_in(2)];

                        % total stage masses (dry+prop)
                        m(:) = [m_in(1) + m_pr(1), m_in(2) + m_pr(2)];

                        % gross initial vehicle mass
                        m_0 = m_o(1);
                end

                function [cost, cost_0] = calculate_stage_costs(obj, m_in)
                        %calculate_stage_costs empirical cost model from inert mass
                        arguments
                                obj TwoStageLV
                                m_in(1, 2) double
                        end
                        cost(:) = [13.52*m_in(1)^(0.55)*1e6, 13.52*m_in(2)^(0.55)*1e6];
                        cost_0 = cost(1) + cost(2);
                end

                function obj = generate_trends(obj, DeltaV, m_pl, delta, g)
                        %generate_trends sweep xs, compute arrays, and track mass/cost minima
                        arguments
                                obj TwoStageLV
                                DeltaV double
                                m_pl double
                                delta(1, 2) double
                                g double
                        end

                        for i=1:100
                                % masses and gross
                                [obj.ms(:, i), obj.m_ins(:, i), obj.m_prs(:, i), obj.m_0s(i)] = ...
                                        obj.calculate_stage_masses(obj.Xs(i), DeltaV, m_pl, delta, g);

                                % costs only when inert masses are physical
                                if (obj.m_ins(1, i) > 1) && (obj.m_ins(2, i) > 1)
                                        [obj.costs(:, i), obj.cost_0s(i)] = obj.calculate_stage_costs(obj.m_ins(:, i));
                                else
                                        obj.costs(:, i) = [nan nan];
                                        obj.cost_0s(i) = nan;
                                end

                                % guard unphysical gross mass
                                if (obj.m_0s(i) < 1)
                                        obj.ms(:, i) = [nan nan];
                                        obj.m_ins(:, i) = [nan nan];
                                        obj.m_prs(:, i) = [nan nan];
                                        obj.m_0s(i) = nan;
                                end

                                % update minima trackers
                                if (obj.m_0s(i) > 1) && (obj.m_0s(i) < obj.m_0_min)
                                        obj.m_0_min = obj.m_0s(i);
                                        obj.X_m_min = obj.Xs(i);
                                        obj.i_m_min = i;
                                end

                                if (obj.cost_0s(i) > 1) && (obj.cost_0s(i) < obj.cost_0_min)
                                        obj.cost_0_min = obj.cost_0s(i);
                                        obj.X_c_min = obj.Xs(i);
                                        obj.i_c_min = i;
                                end
                        end

                        % generate figures (hidden)
                        obj = obj.generate_mass_fig();
                        obj = obj.generate_cost_fig();
                end

                function obj = generate_mass_fig(obj)
                        %generate_mass_fig build hidden figure for mass trends
                        arguments
                                obj TwoStageLV
                        end

                        obj.mass_fig = figure('visible','off');
                        hold on;

                        % masses in tonnes; mark minimum with a dot
                        plot(obj.Xs, obj.m_0s/1e3, Color="red", LineStyle="-", MarkerEdgeColor="black", Marker=".", MarkerIndices=obj.i_m_min);
                        plot(obj.Xs, obj.ms(1, :)/1e3, Color="blue");
                        plot(obj.Xs, obj.ms(2, :)/1e3, Color="green");

                        lgd = legend("Gross Vehicle Mass", "$1^{\mathrm{st}}$ Stage Mass", "$2^{\mathrm{nd}}$ Stage Mass", ...
                                Interpreter="latex", Location="northeastoutside");

                        infoLines = [
                                sprintf("Minimum Gross Vehicle Mass: %.4G $\\mathrm{[t]}$", obj.m_0_min/1e3);
                                sprintf("Min. Mass Program Cost: %.4G $\\mathrm{[\\$B2025]}$", obj.cost_0s(obj.i_m_min)/1e9);
                                sprintf("Min. Mass Stage 1 $\\Delta V$ Fraction: %.4G", obj.X_m_min)
                        ];

                        % layout helpers for the info textbox
                        ax = gca;
                        obj.mass_fig.Units = "normalized";
                        ax.Units  = "normalized";
                        lgd.Units = "normalized";

                        axpos = ax.Position;
                        lgdpos = lgd.Position;
                        gap = 0.010;

                        leftX = axpos(3);
                        bottomY = axpos(2);
                        topY = lgdpos(2) - gap;
                        boxW = max(0.22, 1 - leftX - gap);
                        boxH = max(0.10, topY - bottomY);

                        annotation(obj.mass_fig, textbox=[leftX bottomY boxW boxH], String=infoLines, ...
                                Interpreter="latex", FitBoxToText="on", VerticalAlignment="middle", FontSize=10);

                        xlim([0 1]); xlabel("$X: \Delta V$ Ratio $\left(\frac{\Delta V_1}{\Delta V_2}\right)$", Interpreter="latex");
                        ylim([0 obj.m_0_min/1e3*2]); ylabel("$m \mathrm{[t]}$", Interpreter="latex");

                        set(gca,'TickLabelInterpreter','latex');
                        title(sprintf("Mass Trends for $1^{\\mathrm{st}}$ Stage: %s; $2^{\\mathrm{nd}}$ Stage: %s", ...
                                obj.propellants(1).displayname, obj.propellants(2).displayname), Interpreter="latex");
                        hold off;
                end

                function obj = generate_cost_fig(obj)
                        %generate_cost_fig build hidden figure for cost trends
                        arguments
                                obj TwoStageLV
                        end

                        obj.cost_fig = figure('visible','off');
                        hold on;

                        plot(obj.Xs, obj.cost_0s(:)/1e9, Color="red", LineStyle="-", MarkerEdgeColor="black", Marker=".", MarkerIndices=obj.i_c_min);
                        plot(obj.Xs, obj.costs(1, :)/1e9, Color="blue");
                        plot(obj.Xs, obj.costs(2, :)/1e9, Color="green");

                        lgd = legend("Program Cost", "$1^{\mathrm{st}}$ Stage Cost", "$2^{\mathrm{nd}}$ Stage Cost", ...
                                Interpreter="latex", Location="northeastoutside");

                        infoLines = [
                                sprintf("Minimum Program Cost: %.4G $\\mathrm{[\\$B2025]}$", obj.cost_0_min/1e9);
                                sprintf("Min. Cost Gross Vehicle Mass: %.4G $\\mathrm{[t]}$", obj.m_0s(obj.i_c_min)/1e3);
                                sprintf("Min. Cost Stage 1 $\\Delta V$ Fraction: %.4G", obj.X_c_min)
                        ];

                        ax = gca;
                        obj.cost_fig.Units = "normalized";
                        ax.Units  = "normalized";
                        lgd.Units = "normalized";

                        axpos = ax.Position;
                        lgdpos = lgd.Position;
                        gap = 0.010;

                        leftX = axpos(3);
                        bottomY = axpos(2);
                        topY = lgdpos(2) - gap;
                        boxW = max(0.22, 1 - leftX - gap);
                        boxH = max(0.10, topY - bottomY);

                        annotation(obj.cost_fig, textbox=[leftX bottomY boxW boxH], String=infoLines, ...
                                Interpreter="latex", FitBoxToText="on", VerticalAlignment="middle", FontSize=10);

                        xlim([0 1]); xlabel("$X: \Delta V$ Ratio $\left(\frac{\Delta V_1}{\Delta V_2}\right)$", Interpreter="latex");
                        ylim([0 obj.cost_0_min/1e9*2]); ylabel("$\mathrm{Cost} \mathrm{[\$B2025]}$", Interpreter="latex");

                        title(sprintf("Cost Trends for $1^{\\mathrm{st}}$ Stage: %s; $2^{\\mathrm{nd}}$ Stage: %s", ...
                                obj.propellants(1).displayname, obj.propellants(2).displayname), Interpreter="latex");

                        set(gca,'TickLabelInterpreter','latex');
                        hold off;
                end

                function obj = save_mass_fig(obj)
                        %save_mass_fig persist mass trends figure (fixed size)
                        arguments
                                obj TwoStageLV
                        end
                        filename = sprintf("./output/%2$s/mass/s1 %1$s - s2 %2$s.jpg", obj.propellants(1).name, obj.propellants(2).name);
                        obj.mass_fig.Units = "centimeters";
                        obj.mass_fig.Position = [2 2 24 12];
                        drawnow; saveas(obj.mass_fig, filename);
                end

                function obj = save_cost_fig(obj)
                        %save_cost_fig persist cost trends figure (fixed size)
                        arguments
                                obj TwoStageLV
                        end
                        filename = sprintf("./output/%2$s/cost/s1 %1$s - s2 %2$s.jpg", obj.propellants(1).name, obj.propellants(2).name);
                        obj.cost_fig.Units = "centimeters";
                        obj.cost_fig.Position = [2 2 24 12];
                        drawnow; saveas(obj.cost_fig, filename);
                end
        end
end
