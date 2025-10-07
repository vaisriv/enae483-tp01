% TwoStageLV.m
classdef TwoStageLV
        %TWOSTAGELV two-stage launch vehicle
        % class to compute mass & cost trends for a 2-stage vehicle

        properties
                % propellants: propellant objects for stage1 and stage2 (1x2)
                propellants(1, 2) PropellantMix = PropellantMix("Blank", "Blank", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false)

                % xs: sweep of stage-1 delta-v fraction (1x100)
                Xs(1, 100) double = linspace(0, 1, 100)

                % ms: stage masses (in kg) for stages [1;2] at each xs (2x100)
                ms(2, 100) double = zeros(2, 100)

                % m_ins: inert mass of each stage (2x100)
                m_ins(2, 100) double = zeros(2, 100)

                % m_prs: propellant mass of each stage (2x100)
                m_prs(2, 100) double = zeros(2, 100)

                % m_0s: gross vehicle initial mass for each xs (1x100)
                m_0s(1, 100) double = zeros(1, 100)

                % m_0_min: minimum found gross mass (scalar)
                m_0_min double = inf

                % x_m_min: xs value that produced min gross mass
                X_m_min double = 0

                % i_m_min: index in xs for min gross mass
                i_m_min double = 1

                % costs: stage cost arrays (same indexing as ms) (2x100)
                costs(2, 100) double = zeros(2, 100)

                % cost_0s: total program cost for each xs (1x100)
                cost_0s(1, 100) double = zeros(1, 100)

                % cost_0_min: minimum program cost found (scalar)
                cost_0_min double = inf

                % x_c_min: xs value for minimum cost
                X_c_min double = 0

                % i_c_min: index for minimum cost
                i_c_min double = 1

                % mass_fig: handle to generated mass trends figure
                mass_fig matlab.ui.Figure

                % cost_fig: handle to generated cost trends figure
                cost_fig matlab.ui.Figure

                % subsystem_masses: struct of masses (in kg) for each subsystem
                % for each stage
                stage_subsystem_masses struct = struct( ...
                        "propellant", [0; 0], ...
                        "propellant_tanks", [0; 0], ...
                        "propellant_tank_insulation", [0; 0], ...
                        "engines_or_casing", [0; 0], ...
                        "thrust_structure", [0; 0], ...
                        "gimbals", [0; 0], ...
                        "avionics", [0; 0], ...
                        "wiring", [0; 0], ...
                        "payload_fairing", 0, ...
                        "inter_tank_fairing", [0; 0], ...
                        "inter_stage_fairing", 0, ...
                        "aft_fairing", 0)

                % estimated_m_stages_margin: total of all subsystem masses (in kg) with a
                % 30% margin for each stage
                estimated_m_stages_margin(2, 1) double = [0; 0]

                estimated_m_total_margin double = 0

                estimated_m_stages(2, 1) double = [0; 0]

                estimated_m_total double = 0

                length_total double = 0

                necessary_engine_count double = 0

                estimated_c_stages(2, 1) double = [0; 0]

                estimated_c_total double = 0
        end

        methods
                function obj = TwoStageLV(propellants, DeltaV, m_pl, delta, g)
                        %constructor: store provided propellant mixes (2-element array)
                        arguments
                                propellants(1, 2) PropellantMix
                                DeltaV double
                                m_pl double
                                delta(1, 2) double
                                g double
                        end
                        obj.propellants = propellants;
                        obj = obj.run_part_01(DeltaV, m_pl, delta, g);
                        obj = obj.run_part_02(m_pl, g);
                end

                function obj = run_part_01(obj, DeltaV, m_pl, delta, g)
                        %RUN_PART_01 runs all necessary code for part 01 of
                        %the LV trade study project
                        arguments
                                obj TwoStageLV
                                DeltaV double
                                m_pl double
                                delta(1, 2) double
                                g double
                        end
                        obj = obj.generate_trends(DeltaV, m_pl, delta, g);
                        % obj = obj.save_mass_fig();
                        % obj = obj.save_cost_fig();
                end

                function obj = run_part_02(obj, m_pl, g)
                        %RUN_PART_02 runs all necessary code for part 02 of
                        %the LV trade study project
                        arguments
                                obj TwoStageLV
                                m_pl double
                                g double
                        end

                        obj = obj.calculate_desired_engine_count(m_pl, g);
                        obj = obj.estimate_subsystem_masses(obj.necessary_engine_count, m_pl);
                end

                function [obj, stage_subsystem_masses, estimated_m_stage, estimated_m_nrec, estimated_m_total, l_total] = estimate_subsystem_masses(obj, n, m_pl)
                        arguments
                                obj TwoStageLV
                                n double
                                m_pl double
                        end

                        oxidizer_masses = [
                                (obj.propellants(1).oxidizer_mix_mass_ratio / obj.propellants(1).oxidizer_mix_mass_ratio + obj.propellants(1).fuel_mix_mass_ratio) * obj.m_prs(1, obj.i_m_min);
                                (obj.propellants(2).oxidizer_mix_mass_ratio / obj.propellants(2).oxidizer_mix_mass_ratio + obj.propellants(2).fuel_mix_mass_ratio) * obj.m_prs(2, obj.i_m_min);
                        ];

                        fuel_masses = [
                                (obj.propellants(1).fuel_mix_mass_ratio / obj.propellants(1).oxidizer_mix_mass_ratio + obj.propellants(1).fuel_mix_mass_ratio) * obj.m_prs(1, obj.i_m_min);
                                (obj.propellants(2).fuel_mix_mass_ratio / obj.propellants(2).oxidizer_mix_mass_ratio + obj.propellants(2).fuel_mix_mass_ratio) * obj.m_prs(2, obj.i_m_min);
                        ];

                        oxidizer_volumes = [
                                oxidizer_masses(1)/obj.propellants(1).density_oxidizer;
                                oxidizer_masses(2)/obj.propellants(2).density_oxidizer;
                        ];

                        fuel_volumes = [
                                fuel_masses(1)/obj.propellants(1).density_fuel;
                                fuel_masses(2)/obj.propellants(2).density_fuel;
                        ];

                        r_tank = 9.8; % [m]

                        payload_height = 338/(5*r_tank - 26);
                        payload_area = pi*r_tank*sqrt(r_tank^2+payload_height^2);

                        oxidizer_tank_heights = oxidizer_volumes./(pi*r_tank^2);
                        oxidizer_tank_areas = 2*pi*r_tank*oxidizer_tank_heights;

                        fuel_tank_heights = fuel_volumes./(pi*r_tank^2);
                        fuel_tank_areas = 2*pi*r_tank*fuel_tank_heights;

                        intertank_heights = [
                                1/2*r_tank;
                                1/2*r_tank;
                        ];
                        intertank_areas = 2*pi*r_tank*intertank_heights;

                        interstage_height = 1/2*r_tank;
                        interstage_area = 2*pi*r_tank*interstage_height;

                        aft_height = 1/2*r_tank;
                        aft_area = 2*pi*r_tank*aft_height;

                        ls = oxidizer_tank_heights+fuel_tank_heights+intertank_heights;
                        l_total = payload_height+sum(ls)+interstage_height+aft_height;
                        
                        stage_subsystem_masses.propellant = [
                                obj.m_prs(1, obj.i_m_min);
                                obj.m_prs(2, obj.i_m_min);
                        ];
                        stage_subsystem_masses.propellant_tanks = [
                                obj.propellants(1).propellant_tank_mer_oxidizer*oxidizer_volumes(1) + obj.propellants(1).propellant_tank_mer_fuel*fuel_volumes(1);
                                obj.propellants(2).propellant_tank_mer_oxidizer*oxidizer_volumes(2) + obj.propellants(2).propellant_tank_mer_fuel*fuel_volumes(2);
                        ];
                        stage_subsystem_masses.propellant_tank_insulation = [
                                obj.propellants(1).propellant_tank_insulation_mer_oxidizer*oxidizer_tank_areas(1) + obj.propellants(1).propellant_tank_insulation_mer_fuel*fuel_tank_areas(1);
                                obj.propellants(2).propellant_tank_insulation_mer_oxidizer*oxidizer_tank_areas(2) + obj.propellants(2).propellant_tank_insulation_mer_fuel*fuel_tank_areas(2);
                        ];
                        stage_subsystem_masses.engines_or_casing = [
                                n*calculate_stage_engine_or_casing_mass(obj, 1, obj.propellants(1).thrust_per_motor_stage1, obj.propellants(1).nozzle_expansion_ratio_stage1);
                                n*calculate_stage_engine_or_casing_mass(obj, 2, obj.propellants(2).thrust_per_motor_stage2, obj.propellants(2).nozzle_expansion_ratio_stage2);
                        ];
                        stage_subsystem_masses.thrust_structure = [
                                n*2.55e-4*obj.propellants(1).thrust_per_motor_stage1;
                                n*2.55e-4*obj.propellants(2).thrust_per_motor_stage2;
                        ];
                        stage_subsystem_masses.gimbals = [
                                n*237.8*(obj.propellants(1).thrust_per_motor_stage1/obj.propellants(1).chamber_pressure_stage1)^0.9375;
                                n*237.8*(obj.propellants(2).thrust_per_motor_stage2/obj.propellants(2).chamber_pressure_stage2)^0.9375;
                        ];
                        stage_subsystem_masses.avionics = [
                                10*obj.ms(1, obj.i_m_min)^0.361;
                                10*obj.ms(2, obj.i_m_min)^0.361;
                        ];
                        stage_subsystem_masses.wiring = [
                                1.058*sqrt(obj.ms(1, obj.i_m_min))*ls(1)^(0.25);
                                1.058*sqrt(obj.ms(2, obj.i_m_min))*ls(2)^(0.25);
                        ];
                        stage_subsystem_masses.payload_fairing = 4.95*payload_area^1.15;
                        stage_subsystem_masses.inter_tank_fairing = [
                                4.95*intertank_areas(1)^1.15;
                                4.95*intertank_areas(2)^1.15;
                        ];
                        stage_subsystem_masses.inter_stage_fairing = 4.95*interstage_area^1.15;
                        stage_subsystem_masses.aft_fairing = 4.95*aft_area^1.15;

                        estimated_m_stage = [
                                stage_subsystem_masses.propellant(1)+stage_subsystem_masses.propellant_tanks(1)+stage_subsystem_masses.propellant_tank_insulation(1)+stage_subsystem_masses.engines_or_casing(1)+stage_subsystem_masses.thrust_structure(1)+stage_subsystem_masses.gimbals(1)+stage_subsystem_masses.avionics(1)+stage_subsystem_masses.wiring(1)+stage_subsystem_masses.inter_tank_fairing(1);
                                stage_subsystem_masses.propellant(2)+stage_subsystem_masses.propellant_tanks(2)+stage_subsystem_masses.propellant_tank_insulation(2)+stage_subsystem_masses.engines_or_casing(2)+stage_subsystem_masses.thrust_structure(2)+stage_subsystem_masses.gimbals(2)+stage_subsystem_masses.avionics(2)+stage_subsystem_masses.wiring(2)+stage_subsystem_masses.inter_tank_fairing(2);
                        ];
                        estimated_m_nrec = stage_subsystem_masses.payload_fairing+stage_subsystem_masses.inter_stage_fairing+stage_subsystem_masses.aft_fairing+m_pl;
                        estimated_m_total = sum(estimated_m_stage)+estimated_m_nrec;

                        obj.stage_subsystem_masses = stage_subsystem_masses;
                        obj.estimated_m_stages = estimated_m_stage;
                        obj.estimated_m_total = estimated_m_total;
                        obj.length_total = l_total;
                end

                function m_e_or_c = calculate_stage_engine_or_casing_mass(obj, stage, T_N, NER)
                        arguments
                                obj TwoStageLV
                                stage double
                                T_N double
                                NER double
                        end
                        if obj.propellants(stage).is_solid
                                m_e_or_c = 0.135*obj.stage_subsystem_masses.propellant(stage);
                        else
                                m_e_or_c = 7.81e-4*T_N + 3.37e-5*T_N*sqrt(NER) + 59;
                        end                 
                end

                function obj = calculate_desired_engine_count(obj, m_pl, g)
                        arguments
                                obj TwoStageLV
                                m_pl double
                                g double
                        end
                        n = 1;
                        valid_engine_count = false;
                        while not(valid_engine_count)
                                [~, ~, est_m_stage, est_m_nrec, est_m_total, ~] = obj.estimate_subsystem_masses(n, m_pl);
                                T_des = [
                                        est_m_total*g/(1.3*n);
                                        (est_m_stage(2)+est_m_nrec)*g/(0.76*n);   
                                ];

                                if (obj.propellants(1).thrust_per_motor_stage1 >= T_des(1)) && (obj.propellants(2).thrust_per_motor_stage2 >= T_des(2))
                                        valid_engine_count = true;
                                else
                                        n = n+1;
                                end
                        end

                        obj.necessary_engine_count = n;
                end

                function obj = calculate_margin_masses_and_costs(obj)
                        arguments
                                obj TwoStageLV
                        end

                        obj.estimated_m_total_margin = obj.estimated_m_total*1.3;
                        obj.estimated_m_stages_margin = [
                                obj.estimated_m_stages(1)*1.3;
                                obj.estimated_m_stages(2)*1.3;
                        ];

                        obj.estimated_c_stages = [
                                13.52*obj.estimated_m_stages_margin(1)^(0.55)*1e6;
                                13.52*obj.estimated_m_stages_margin(2)^(0.55)*1e6;
                        ];
                        obj.estimated_c_total = 13.52*obj.estimated_m_total_margin^(0.55)*1e6;
                end

                function [m, m_in, m_pr, m_0] = calculate_stage_masses(obj, X, DeltaV, m_pl, delta, g)
                        %CALCULATE_STAGE_MASSES compute masses for a given split X
                        %
                        % inputs:
                        %   x       - fraction of total delta-v assigned to stage1
                        %   deltav  - total mission delta-v (same units as ve*ln)
                        %   m_pl    - payload mass (scalar)
                        %   delta   - inert mass fraction vector for stages [1,2]
                        %   g       - gravitational accel used to convert isp->ve
                        %
                        % outputs:
                        %   m       - stage dry+prop masses [m1; m2]
                        %   m_in    - inert masses for each stage
                        %   m_pr    - propellant masses for each stage
                        %   m_0     - gross initial vehicle mass (stage1 initial mass)
                        arguments
                                obj TwoStageLV
                                X double
                                DeltaV double
                                m_pl double
                                delta(1, 2) double
                                g double
                        end

                        % effective exhaust velocities for each stage: ve = isp * g
                        Ve(:) = [
                                g*obj.propellants(1).specific_impulse ...
                                g*obj.propellants(2).specific_impulse
                        ];

                        % split the total delta-v into stage1 and stage2 pieces
                        DV(:) = [
                                X*DeltaV ...
                                (1-X)*DeltaV
                        ];

                        % mass-ratio exponentials for each stage: r = exp(-DV/Ve)
                        r(:) = [
                                exp(-DV(1)/Ve(1)) ...
                                exp(-DV(2)/Ve(2))
                        ];

                        % lambda = r - inert_fraction (used to algebraically solve masses)
                        lambda(:) = [
                                r(1) - delta(1) ...
                                r(2) - delta(2)
                        ];

                        % solve upward from payload: stage2 initial mass m_o2 = m_pl / lambda2
                        m_o(2) = m_pl/lambda(2);
                        % stage1 initial mass found by dividing by lambda1
                        m_o(1) = m_o(2)/lambda(1);

                        % inert (structural) mass for each stage = delta * m_o(stage)
                        m_in(:) = [
                                delta(1)*m_o(1) ...
                                delta(2)*m_o(2)
                        ];

                        % propellant masses: stage1 prop = m_o1 - m_o2 - m_in1
                        m_pr(:) = [
                                m_o(1) - m_o(2) - m_in(1) ...
                                m_o(2) - m_pl - m_in(2)
                        ];

                        % total stage mass = inert + propellant
                        m(:) = [
                                m_in(1) + m_pr(1) ...
                                m_in(2) + m_pr(2)
                        ];

                        % overall initial mass is stage1 initial mass
                        m_0 = m_o(1);
                end

                function [cost, cost_0] = calculate_stage_costs(obj, m_in)
                        %CALCULATE_STAGE_COSTS empirical cost model based on inert mass
                        %
                        % inputs:
                        %   m_in - inert mass vector [m_in1; m_in2]
                        %
                        % outputs:
                        %   cost   - cost per stage (same units as formula)
                        %   cost_0 - total program cost (sum)
                        arguments
                                obj TwoStageLV
                                m_in(1, 2) double
                        end

                        % simple power-law cost scaling, multiplied by 1e6 to set units
                        cost(:) = [
                                13.52*m_in(1)^(0.55)*1e6 ...
                                13.52*m_in(2)^(0.55)*1e6
                        ];

                        % total program cost
                        cost_0 = cost(1) + cost(2);
                end

                function obj = generate_trends(obj, DeltaV, m_pl, delta, g)
                        %GENERATE_TRENDS sweep xs, compute mass & cost trends, and find minima
                        %
                        % performs 100 evals over obj.xs, populates object arrays and figs
                        arguments
                                obj TwoStageLV
                                DeltaV double
                                m_pl double
                                delta(1, 2) double
                                g double
                        end

                        for i=1:100
                                % compute masses for this xs
                                [obj.ms(:, i), obj.m_ins(:, i), obj.m_prs(:, i), obj.m_0s(i)] = obj.calculate_stage_masses(obj.Xs(i), DeltaV, m_pl, delta, g);

                                % only compute cost if inert masses look reasonable (>1 kg)
                                if (obj.m_ins(1, i) > 1) && (obj.m_ins(2, i) > 1)
                                        [obj.costs(:, i), obj.cost_0s(i)] = obj.calculate_stage_costs(obj.m_ins(:, i));
                                else
                                        % otherwise mark cost as nan
                                        obj.costs(:, i) = [nan nan];
                                        obj.cost_0s(i) = nan;
                                end

                                % if gross mass is unphysical (<1) mark all masses nan
                                if (obj.m_0s(i) < 1)
                                        obj.ms(:, i) = [nan nan];
                                        obj.m_ins(:, i) = [nan nan];
                                        obj.m_prs(:, i) = [nan nan];
                                        obj.m_0s(i) = nan;
                                end

                                % update minimum gross mass tracker
                                if (obj.m_0s(i) > 1) && (obj.m_0s(i) < obj.m_0_min)
                                        obj.m_0_min = obj.m_0s(i);
                                        obj.X_m_min = obj.Xs(i);
                                        obj.i_m_min = i;
                                end

                                % update minimum cost tracker
                                if (obj.cost_0s(i) > 1) && (obj.cost_0s(i) < obj.cost_0_min)
                                        obj.cost_0_min = obj.cost_0s(i);
                                        obj.X_c_min = obj.Xs(i);
                                        obj.i_c_min = i;
                                end
                        end

                        % generate figure objects (saved in obj)
                        obj = obj.generate_mass_fig();
                        obj = obj.generate_cost_fig();
                end

                function obj = generate_mass_fig(obj)
                        %GENERATE_MASS_FIG create a hidden figure summarizing mass trends
                        arguments
                                obj TwoStageLV
                        end

                        obj.mass_fig = figure('visible','off');
                        hold on;

                        % plot gross vehicle mass and per-stage masses (converted to tonnes)
                        plot(obj.Xs, obj.m_0s/1e3, Color="red", LineStyle="-", MarkerEdgeColor="black", Marker=".", MarkerIndices=obj.i_m_min);
                        plot(obj.Xs, obj.ms(1, :)/1e3, Color="blue");
                        plot(obj.Xs, obj.ms(2, :)/1e3, Color="green");

                        % legend with latex interpreter
                        lgd = legend( ...
                                "Gross Vehicle Mass", ...
                                "$1^{\mathrm{st}}$ Stage Mass", ...
                                "$2^{\mathrm{nd}}$ Stage Mass", ...
                                Interpreter="latex", ...
                                Location="northeastoutside");

                        % info lines displayed in an annotation textbox (uses latex formatting)
                        infoLines = [
                                sprintf("Minimum Gross Vehicle Mass: %.4G $\\mathrm{[t]}$", obj.m_0_min/1e3);
                                sprintf("Min. Mass Program Cost: %.4G $\\mathrm{[\\$B2025]}$", obj.cost_0s(obj.i_m_min)/1e9);
                                sprintf("Min. Mass Stage 1 $\\Delta V$ Fraction: %.4G", obj.X_m_min)
                        ];

                        % normalize positions so annotation can be placed outside the axis
                        ax = gca;
                        obj.mass_fig.Units = "normalized";
                        ax.Units  = "normalized";
                        lgd.Units = "normalized";

                        % compute textbox position to the right of the axes and below legend
                        axpos = ax.Position;
                        lgdpos = lgd.Position;
                        gap = 0.010;

                        leftX = axpos(3); % x start just to the right of axis width
                        bottomY = axpos(2); % same bottom as axis
                        topY = lgdpos(2) - gap; % top limited by legend top minus gap
                        boxW = max(0.22, 1 - leftX - gap); % ensure textbox width reasonable
                        boxH = max(0.10, topY - bottomY); % ensure textbox height reasonable

                        % add annotation textbox with the info lines
                        annotation(obj.mass_fig, ...
                                textbox=[leftX bottomY boxW boxH], ...
                                String=infoLines, ...
                                Interpreter="latex", ...
                                FitBoxToText="on", ...
                                VerticalAlignment="middle", ...
                                FontSize=10);

                        % axis labels and limits
                        xlim([0 1]);
                        xlabel("$X$", Interpreter="latex");

                        ylim([0 obj.m_0_min/1e3*2]); % y-lim set to twice the minimum gross mass (tonnes)
                        ylabel("$m \mathrm{[t]}$", Interpreter="latex");

                        % set title using propellant display names
                        set(gca,'TickLabelInterpreter','latex');
                        title(sprintf( ...
                                "Mass Trends for $1^{\\mathrm{st}}$ Stage: %s; $2^{\\mathrm{nd}}$ Stage: %s", ...
                                obj.propellants(1).displayname, ...
                                obj.propellants(2).displayname), ...
                                Interpreter="latex");

                        hold off;
                end

                function obj = generate_cost_fig(obj)
                        %GENERATE_COST_FIG create hidden figure summarizing cost trends
                        arguments
                                obj TwoStageLV
                        end

                        obj.cost_fig = figure('visible','off');
                        hold on;

                        % plot total program cost and per-stage costs (converted to billion dollars)
                        plot(obj.Xs, obj.cost_0s(:)/1e9, Color="red", LineStyle="-", MarkerEdgeColor="black", Marker=".", MarkerIndices=obj.i_c_min);
                        plot(obj.Xs, obj.costs(1, :)/1e9, Color="blue");
                        plot(obj.Xs, obj.costs(2, :)/1e9, Color="green");

                        % legend and info textbox similar to mass figure
                        lgd = legend( ...
                                "Program Cost", ...
                                "$1^{\mathrm{st}}$ Stage Cost", ...
                                "$2^{\mathrm{nd}}$ Stage Cost", ...
                                Interpreter="latex", ...
                                Location="northeastoutside");

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

                        annotation(obj.cost_fig, ...
                                textbox=[leftX bottomY boxW boxH], ...
                                String=infoLines, ...
                                Interpreter="latex", ...
                                FitBoxToText="on", ...
                                VerticalAlignment="middle", ...
                                FontSize=10);

                        xlim([0 1]);
                        xlabel("$X$", Interpreter="latex")
                        ylim([0 obj.cost_0_min/1e9*2]); % ylim set to twice the minimum program cost (billion $)
                        ylabel("$\mathrm{Cost} \mathrm{[\$B2025]}$", Interpreter="latex")

                        title(sprintf( ...
                                "Cost Trends for $1^{\\mathrm{st}}$ Stage: %s; $2^{\\mathrm{nd}}$ Stage: %s", ...
                                obj.propellants(1).displayname, ...
                                obj.propellants(2).displayname), ...
                                Interpreter="latex")

                        set(gca,'TickLabelInterpreter','latex');
                        hold off;
                end

                function obj = save_mass_fig(obj)
                        %SAVE_MASS_FIG save mass figure to disk with fixed dimensions
                        arguments
                                obj TwoStageLV
                        end

                        % filename uses propellant short names (stage1, stage2) in folder structure
                        filename = sprintf("./images/mass/%2$s/s1 %1$s - s2 %2$s.jpg", obj.propellants(1).name, obj.propellants(2).name);

                        % set figure size in centimeters then save
                        obj.mass_fig.Units = "centimeters";
                        obj.mass_fig.Position = [2 2 24 12];

                        drawnow; % ensure rendering is flushed
                        saveas(obj.mass_fig, filename);
                end

                function obj = save_cost_fig(obj)
                        %SAVE_COST_FIG save cost figure to disk with fixed dimensions
                        arguments
                                obj TwoStageLV
                        end

                        filename = sprintf("./images/cost/%2$s/s1 %1$s - s2 %2$s.jpg", obj.propellants(1).name, obj.propellants(2).name);

                        obj.cost_fig.Units = "centimeters";
                        obj.cost_fig.Position = [2 2 24 12];

                        drawnow;
                        saveas(obj.cost_fig, filename);
                end
        end
end