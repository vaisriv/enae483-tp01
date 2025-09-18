classdef TwoStageLV
        %TWOSTAGELV Two-Stage Launch Vehicle
        % See also PROPELLANTMIX

        properties
                % PROPELLANTS Propellants used in Stages 1 and 2
                propellants(1, 2) PropellantMix = PropellantMix("Blank", "Blank", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

                % XS Array of all X values
                Xs(1, 100) double = linspace(0, 1, 100)

                % MS Array of mass values
                ms(2, 100) double = zeros(2, 100)

                % M_INS Array of inert mass values
                m_ins(2, 100) double = zeros(2, 100)

                % M_PRS Array of propellant mass values
                m_prs(2, 100) double = zeros(2, 100)

                % M_0S Array of total mass values
                m_0s(1, 100) double = zeros(1, 100)

                % M_0_MIN Minimum Total Mass
                m_0_min double = inf

                % X_M_MIN DeltaV Fraction for Minimum Total Mass
                X_m_min double = 0

                % I_M_MIN Index for Minimum Total Mass
                i_m_min double = 1

                % COSTS Array of cost values
                costs(2, 100) double = zeros(2, 100)

                % COST_0S Array of total cost values
                cost_0s(1, 100) double = zeros(1, 100)

                % COST_0_MIN Minimum Total Cost
                cost_0_min double = inf

                % X_C_MIN DeltaV Fraction for Minimum Total Cost
                X_c_min double = 0

                % I_C_MIN Index for Minimum Total Cost
                i_c_min double = 1

                % MASS_FIG Mass Trends Figure
                mass_fig matlab.ui.Figure

                % COST_FIG Cost Trends Figure
                cost_fig matlab.ui.Figure
        end

        methods
                function obj = TwoStageLV(propellants)
                        %TWOSTAGELV Construct an instance of this class
                        arguments
                                % PROPELLANTS Propellants used in Stages 1 and 2
                                propellants(1, 2) PropellantMix
                        end
                        obj.propellants = propellants;
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
                                g double = 9.81
                        end

                        Ve(:) = [
                                g*obj.propellants(1).specific_impulse ...
                                g*obj.propellants(2).specific_impulse
                                ];

                        DV(:) = [
                                X*DeltaV ...
                                (1-X)*DeltaV
                                ];

                        r(:) = [
                                exp(-DV(1)/Ve(1)) ...
                                exp(-DV(2)/Ve(2))
                                ];

                        lambda(:) = [
                                r(1) - delta(1) ...
                                r(2) - delta(2)
                                ];

                        m_o(2) = m_pl/lambda(2);
                        m_o(1) = m_o(2)/lambda(1);

                        m_in(:) = [
                                delta(1)*m_o(1) ...
                                delta(2)*m_o(2)
                                ];

                        m_pr(:) = [
                                m_o(1) - m_o(2) - m_in(1) ...
                                m_o(2) - m_pl - m_in(2)
                                ];

                        m(:) = [
                                m_in(1) + m_pr(1) ...
                                m_in(2) + m_pr(2)
                                ];

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

                        cost(:) = [
                                13.52*m_in(1)^(0.55)*1e6 ...
                                13.52*m_in(2)^(0.55)*1e6
                                ];

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
                                g double = 9.81
                        end

                        for i=1:100
                                [obj.ms(:, i), obj.m_ins(:, i), obj.m_prs(:, i), obj.m_0s(i)] = obj.calculate_stage_masses(obj.Xs(i), DeltaV, m_pl, delta, g);

                                if (obj.m_ins(1, i) > 1) && (obj.m_ins(2, i) > 1)
                                        [obj.costs(:, i), obj.cost_0s(i)] = obj.calculate_stage_costs(obj.m_ins(:, i));
                                else
                                        obj.costs(:, i) = [nan nan];
                                        obj.cost_0s(i) = nan;
                                end

                                if (obj.m_0s(i) < 1)
                                        obj.ms(:, i) = [nan nan];
                                        obj.m_ins(:, i) = [nan nan];
                                        obj.m_prs(:, i) = [nan nan];
                                        obj.m_0s(i) = nan;
                                end

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

                        obj = obj.generate_mass_fig();
                        obj = obj.generate_cost_fig();
                end

                function obj = generate_mass_fig(obj)
                        %GENERATE_MASS_FIG Generate mass figure
                        arguments
                                % OBJ The object itself
                                obj TwoStageLV
                        end

                        obj.mass_fig = figure('visible','off');
                        hold on;

                        plot(obj.Xs, obj.m_0s/1e3, Color="red", LineStyle="-", MarkerEdgeColor="black", Marker=".", MarkerIndices=obj.i_m_min);
                        plot(obj.Xs, obj.ms(1, :)/1e3, Color="blue");
                        plot(obj.Xs, obj.ms(2, :)/1e3, Color="green");

                        lgd = legend( ...
                                "Gross Vehicle Mass", ...
                                "$1^{\mathrm{st}}$ Stage Mass", ...
                                "$2^{\mathrm{nd}}$ Stage Mass", ...
                                Interpreter="latex", ...
                                Location="northeastoutside");

                        infoLines = [
                                sprintf("Minimum Gross Vehicle Mass: %.4G $\\mathrm{[t]}$", obj.m_0_min/1e3);
                                sprintf("Min. Mass Program Cost: %.4G $\\mathrm{[\\$B2025]}$", obj.cost_0s(obj.i_m_min)/1e9);
                                sprintf("Min. Mass Stage 1 $\\Delta V$ Fraction: %.4G", obj.X_m_min)
                                ];

                        ax = gca;
                        obj.mass_fig.Units = "normalized";
                        ax.Units  = "normalized";
                        lgd.Units = "normalized";

                        axpos = ax.Position;
                        lgdpos = lgd.Position;
                        gap   = 0.010;

                        leftX   = max(axpos(1) + axpos(3) + gap, lgdpos(1));
                        bottomY = axpos(2);
                        topY    = lgdpos(2) - gap;
                        boxW    = max(0.22, 1 - leftX - gap);
                        boxH    = max(0.10, topY - bottomY);

                        ann = annotation(obj.mass_fig, ...
                                textbox=[leftX bottomY boxW boxH], ...
                                String=infoLines, ...
                                Interpreter="latex", ...
                                FitBoxToText="on", ...
                                VerticalAlignment="middle", ...
                                FontSize=10);

                        xlim([0 1]);
                        xlabel("$X$", Interpreter="latex");

                        ylim([0 obj.m_0_min/1e3*2]);
                        ylabel("$m \mathrm{[t]}$", Interpreter="latex");

                        set(gca,'TickLabelInterpreter','latex');
                        title(sprintf( ...
                                "Mass Trends for $1^{\\mathrm{st}}$ Stage: %s; $2^{\\mathrm{nd}}$ Stage: %s", ...
                                obj.propellants(1).displayname, ...
                                obj.propellants(2).displayname), ...
                                Interpreter="latex");

                        hold off;
                end

                function obj = generate_cost_fig(obj)
                        %GENERATE_COST_FIG Generate cost figure
                        arguments
                                % OBJ The object itself
                                obj TwoStageLV
                        end

                        obj.cost_fig = figure('visible','off');
                        hold on;

                        plot(obj.Xs, obj.cost_0s(:)/1e9, Color="red", LineStyle="-", MarkerEdgeColor="black", Marker=".", MarkerIndices=obj.i_c_min);
                        plot(obj.Xs, obj.costs(1, :)/1e9, Color="blue");
                        plot(obj.Xs, obj.costs(2, :)/1e9, Color="green");

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

                        leftX   = max(axpos(1) + axpos(3) + gap, lgdpos(1));
                        bottomY = axpos(2);
                        topY    = lgdpos(2) - gap;
                        boxW    = max(0.22, 1 - leftX - gap);
                        boxH    = max(0.10, topY - bottomY);

                        ann = annotation(obj.cost_fig, ...
                                textbox=[leftX bottomY boxW boxH], ...
                                String=infoLines, ...
                                Interpreter="latex", ...
                                FitBoxToText="on", ...
                                VerticalAlignment="middle", ...
                                FontSize=10);

                        xlim([0 1]);
                        xlabel("$X$", Interpreter="latex")
                        ylim([0 obj.cost_0_min/1e9*2]);
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
                        %SAVE_MASS_FIG Save mass figure
                        arguments
                                % OBJ The object itself
                                obj TwoStageLV
                        end

                        filename = sprintf("./images/mass/%2$s/s1 %1$s - s2 %2$s.jpg", obj.propellants(1).name, obj.propellants(2).name);

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

                        filename = sprintf("./images/cost/%2$s/s1 %1$s - s2 %2$s.jpg", obj.propellants(1).name, obj.propellants(2).name);

                        obj.cost_fig.Units = "centimeters";
                        obj.cost_fig.Position = [2 2 32 24];

                        drawnow;
                        saveas(obj.cost_fig, filename);
                end
        end
end
