classdef TwoStageLV
        %TWOSTAGELV Two-Stage Launch Vehicle
        % See also PROPELLANTMIX

        properties
                % PROPELLANTS Propellants used in Stages 1 and 2
                propellants(1, 2) PropellantMix = PropellantMix("Blank", "Blank", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

                % M_0_MIN Minimum Total Mass
                m_0_min double

                % X_M_MIN DeltaV Fraction for Minimum Total Mass
                X_m_min double

                % I_M_MIN Index for Minimum Total Mass
                i_m_min double

                % COST_0_MIN Minimum Total Cost
                cost_0_min double

                % X_C_MIN DeltaV Fraction for Minimum Total Cost
                X_c_min double

                % I_C_MIN Index for Minimum Total Cost
                i_c_min double

                % MASS_FIG Mass Trends Figure
                mass_fig matlab.ui.Figure

                % COST_FIG Cost Trends Figure
                cost_fig matlab.ui.Figure
        end

        methods
                function obj = TwoStageLV(stage1_propellant, stage2_propellant)
                        %TWOSTAGELV Construct an instance of this class
                        arguments
                                % STAGE1_PROPELLANT Propellant used in Stage 1
                                stage1_propellant PropellantMix

                                % STAGE2_PROPELLANT Propellant used in Stage 2
                                stage2_propellant PropellantMix
                        end
                        obj.propellants(1) = stage1_propellant;
                        obj.propellants(2) = stage2_propellant;
                end

                function [m, m_in, m_pr, m_0] = calculate_stage_masses(obj, X, DeltaV, m_pl, delta, g)
                        %CALCULATE_STAGE_MASSES Calculate stage masses
                        arguments
                                % OBJ The object itself
                                obj TwoStageLV

                                % X Fraction of Delta V
                                X double

                                % DELTAV Desired Delta V
                                DeltaV double

                                % M_PL Mass of Payload
                                m_pl double

                                % DELTA delta
                                delta(1, 2) double

                                % G Acceleration due to gravity on Earth
                                g double = 9.82
                        end

                        Ve(1) = g*obj.propellants(1).specific_impulse;
                        Ve(2) = g*obj.propellants(2).specific_impulse;

                        DV(1) = X*DeltaV;
                        DV(2) = (1-X)*DeltaV;

                        r(1) = exp(-DV(1)/Ve(1));
                        r(2) = exp(-DV(2)/Ve(2));

                        lambda(1) = r(1) - delta(1);
                        lambda(2) = r(2) - delta(2);

                        m_o(2) = m_pl/lambda(2);
                        m_o(1) = m_o(2)/lambda(1);

                        m_in(1) = delta(1)*m_o(1);
                        m_in(2) = delta(2)*m_o(2);

                        m_pr(1) = m_o(1) - m_o(2) - m_in(1);
                        m_pr(2) = m_o(2) - m_pl - m_in(2);

                        m(1) = m_in(1) + m_pr(1);
                        m(2) = m_in(2) + m_pr(2);
                        m_0 = m_o(1);
                end

                function [cost, cost_0] = calculate_stage_costs(obj, m_in)
                        %CALCULATE_STAGE_COSTS Calculate stage costs
                        arguments
                                % OBJ The object itself
                                obj TwoStageLV

                                % M_IN Inert Mass of Stages
                                m_in(1, 2) double
                        end

                        cost(1) = 13.52*m_in(1)^(0.55)*1e6;
                        cost(2) = 13.52*m_in(2)^(0.55)*1e6;

                        cost_0 = cost(1) + cost(2);
                end

                function obj = generate_trends(obj, DeltaV, m_pl, delta, g)
                        %GENERATE_TRENDS Generate and plot trends
                        arguments
                                % OBJ The object itself
                                obj TwoStageLV

                                % DELTAV Desired Delta V
                                DeltaV double

                                % M_PL Mass of Payload
                                m_pl double

                                % DELTA delta
                                delta(1, 2) double

                                % G Acceleration due to gravity on Earth
                                g double = 9.82
                        end

                        Xs = linspace(0, 1, 100);
                        ms = zeros(2, 100);
                        m_ins = zeros(2, 100);
                        m_prs = zeros(2, 100);
                        m_0s = zeros(1, 100);

                        costs = zeros(2, 100);
                        cost_0s = zeros(1, 100);

                        obj.m_0_min = inf;
                        obj.i_m_min = 0;

                        obj.cost_0_min = inf;
                        obj.i_c_min = 0;

                        for i=1:100
                                [ms(:, i), m_ins(:, i), m_prs(:, i), m_0s(i)] = obj.calculate_stage_masses(Xs(i), DeltaV, m_pl, delta, g);

                                if (m_ins(1, i) > 1) && (m_ins(2, i) > 1)
                                        [costs(:, i), cost_0s(i)] = obj.calculate_stage_costs(m_ins(:, i));
                                else
                                        costs(1, i) = nan;
                                        costs(2, i) = nan;
                                        cost_0s(i) = nan;
                                end

                                if (m_0s(i) < 1)
                                        ms(1, i) = nan;
                                        ms(2, i) = nan;
                                        m_ins(1, i) = nan;
                                        m_ins(2, i) = nan;
                                        m_prs(1, i) = nan;
                                        m_prs(2, i) = nan;
                                        m_0s(i) = nan;
                                end

                                if (m_0s(i) > 1) && (m_0s(i) < obj.m_0_min)
                                        obj.m_0_min = m_0s(i);
                                        obj.X_m_min = Xs(i);
                                        obj.i_m_min = i;
                                end

                                if (cost_0s(i) > 1) && (cost_0s(i) < obj.cost_0_min)
                                        obj.cost_0_min = cost_0s(i);
                                        obj.X_c_min = Xs(i);
                                        obj.i_c_min = i;
                                end
                        end

                        obj = obj.generate_mass_fig(Xs, ms, m_0s);
                        obj = obj.generate_cost_fig(Xs, costs, cost_0s);
                end

                function obj = generate_mass_fig(obj, Xs, ms, m_0s)
                        %GENERATE_MASS_FIG Generate mass figure
                        arguments
                                % OBJ The object itself
                                obj TwoStageLV

                                % XS Array of X values
                                Xs(1, 100) double

                                % MS Array of mass values
                                ms(2, 100) double

                                % M_0S Array of total mass values
                                m_0s(1, 100) double
                        end

                        obj.mass_fig = figure('visible','off');
                        hold on;

                        % plot(Xs, m_0s/1000, Color="red", LineStyle="-", MarkerEdgeColor="black", Marker=".", MarkerIndices=obj.i_m_min);

                        plot(Xs, m_0s/1000, Color="red");
                        plot(Xs, ms(1, :)/1000, Color="blue");
                        plot(Xs, ms(2, :)/1000, Color="green");

                        plot(obj.X_m_min, obj.m_0_min/1000, Color="none", MarkerEdgeColor="black", Marker=".");

                        legend( ...
                                "Gross Vehicle Mass", ...
                                "$1^{\mathrm{st}}$ Stage Mass", ...
                                "$2^{\mathrm{nd}}$ Stage Mass", ...
                                sprintf("Minimum Gross Vehicle Mass: %.4G $\\mathrm{[t]}$", obj.m_0_min/1000), ...
                                Interpreter="latex", ...
                                Location="eastoutside");

                        xlim([0 1]);
                        xlabel("$X$", Interpreter="latex")

                        ylim([0 obj.m_0_min/1000*2]);
                        ylabel("$m \mathrm{[t]}$", Interpreter="latex")

                        title(sprintf( ...
                                        "Mass Trends for $1^{\\mathrm{st}}$ Stage: %s; $2^{\\mathrm{nd}}$ Stage: %s", ...
                                        obj.propellants(1).displayname, ...
                                        obj.propellants(2).displayname), ...
                                Interpreter="latex")

                        hold off;
                end

                function obj = generate_cost_fig(obj, Xs, costs, cost_0s)
                        %GENERATE_COST_FIG Generate cost figure
                        arguments
                                % OBJ The object itself
                                obj TwoStageLV

                                % XS Array of X values
                                Xs(1, 100) double

                                % COSTS Array of cost values
                                costs(2, 100) double

                                % COST_0S Array of total cost values
                                cost_0s(1, 100) double
                        end

                        obj.cost_fig = figure('visible','off');
                        hold on;

                        plot(Xs, cost_0s(:)/1e9, Color="red");
                        plot(Xs, costs(1, :)/1e9, Color="blue");
                        plot(Xs, costs(2, :)/1e9, Color="green");

                        plot(obj.X_c_min, obj.cost_0_min/1e9, Color="none", MarkerEdgeColor="black", Marker=".");

                        legend( ...
                                "Gross Vehicle Cost", ...
                                "$1^{\mathrm{st}}$ Stage Cost", ...
                                "$2^{\mathrm{nd}}$ Stage Cost", ...
                                sprintf("Minimum Gross Vehicle Cost: %.4G $\\mathrm{[\\$B2025]}$", obj.cost_0_min/1e9), ...
                                Interpreter="latex", ...
                                Location="eastoutside");

                        xlim([0 1]);
                        xlabel("$X$", Interpreter="latex")

                        ylim([0 obj.cost_0_min/1e9*2]);
                        ylabel("$\mathrm{Cost} \mathrm{[\$B2025]}$", Interpreter="latex")

                        title(sprintf( ...
                                        "Cost Trends for $1^{\\mathrm{st}}$ Stage: %s; $2^{\\mathrm{nd}}$ Stage: %s", ...
                                        obj.propellants(1).displayname, ...
                                        obj.propellants(2).displayname), ...
                                Interpreter="latex")

                        hold off;
                end

                function obj = save_mass_fig(obj)
                        %SAVE_MASS_FIG Save mass figure
                        arguments
                                % OBJ The object itself
                                obj TwoStageLV
                        end

                        filename = sprintf("./images/mass/s1 %s - s2 %s.jpg", obj.propellants(1).name, obj.propellants(2).name);

                        obj.mass_fig.Units = "centimeters";
                        obj.mass_fig.Position = [2 2 32 24];

                        drawnow;
                        saveas(obj.mass_fig, filename);
                end

                function obj = save_cost_fig(obj)
                        %SAVE_COST_FIG Save cost figure
                        arguments
                                % OBJ The object itself
                                obj TwoStageLV
                        end

                        filename = sprintf("./images/cost/s1 %s - s2 %s.jpg", obj.propellants(1).name, obj.propellants(2).name);

                        obj.cost_fig.Units = "centimeters";
                        obj.cost_fig.Position = [2 2 32 24];

                        drawnow;
                        saveas(obj.cost_fig, filename);
                end
        end
end
