% TwoStageLV.m
classdef TwoStageLV
        %TWOSTAGELV two-stage launch vehicle
        % class to compute mass & cost trends for a 2-stage vehicle

        properties
                % propellants: propellant objects for stage1 and stage2 (1x2)
                propellants(1, 2) PropellantMix = PropellantMix("Blank", "Blank", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

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
        end

        methods
                function obj = TwoStageLV(propellants)
                        %constructor: store provided propellant mixes (2-element array)
                        arguments
                                propellants(1, 2) PropellantMix
                        end
                        obj.propellants = propellants;
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
                        gap   = 0.010;

                        leftX   = axpos(3);               % x start just to the right of axis width
                        bottomY = axpos(2);               % same bottom as axis
                        topY    = lgdpos(2) - gap;        % top limited by legend top minus gap
                        boxW    = max(0.22, 1 - leftX - gap); % ensure textbox width reasonable
                        boxH    = max(0.10, topY - bottomY);  % ensure textbox height reasonable

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
                        gap   = 0.010;

                        leftX   = axpos(3);
                        bottomY = axpos(2);
                        topY    = lgdpos(2) - gap;
                        boxW    = max(0.22, 1 - leftX - gap);
                        boxH    = max(0.10, topY - bottomY);

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

                        drawnow;           % ensure rendering is flushed
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