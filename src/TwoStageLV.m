classdef TwoStageLV
        %TWOSTAGELV Two-Stage Launch Vehicle
        % See also PROPELLANTMIX

        properties
                % PROPELLANTS Propellants used in Stages 1 and 2
                propellants(1, 2) PropellantMix = PropellantMix("Blank", "Blank", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

                % M_0_MIN Minimum Total Mass
                m_0_min double

                % B2025_MIN Minimum Cost
                B2025_min double

                % FIG Mass Trends Figure
                mass_fig matlab.ui.Figure
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

                function [m, m_pr, m_0] = calculate_stage_masses(obj, X, DeltaV, m_pl, delta, g)
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
                        Ve(2) = g*obj.propellants(1).specific_impulse;

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

                function obj = generate_mass_trends(obj, DeltaV, m_pl, delta, g)
                        %GENERATE_MASS_TRENDS Generate and plot mass trends
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

                        % setup
                        Xs = linspace(0, 1, 100);
                        ms = zeros(2, 100);
                        m_prs = zeros(2, 100);
                        m_0s = zeros(1, 100);

                        obj.m_0_min = inf;
                        min_index = 0;

                        % calculating masses for each fraction
                        for i=1:100
                                [ms(:, i), m_prs(:, i), m_0s(i)] = obj.calculate_stage_masses(Xs(i), DeltaV, m_pl, delta, g);
                                if m_0s(i)>1 && m_0s(i)<obj.m_0_min
                                        obj.m_0_min = m_0s(i);
                                        min_index = i;
                                end
                        end

                        % figure
                        obj.mass_fig = figure;
                        hold on;

                        title(sprintf( ...
                                        "$1^{\\mathrm{st}}$ Stage: %s; $2^{\\mathrm{nd}}$ Stage: %s", ...
                                        obj.propellants(1).displayname, ...
                                        obj.propellants(2).displayname), ...
                                Interpreter="latex")
                        xlabel("$X$", Interpreter="latex")
                        ylabel("$m \mathrm{[t]}$", Interpreter="latex")

                        plot(Xs, ms(1, :)/1000, Color="red");
                        plot(Xs, ms(2, :)/1000, Color="green");
                        plot(Xs, m_0s/1000, Color="blue");
                        plot(Xs(min_index), obj.m_0_min/1000, Color="black", Marker="o");
                        legend( ...
                                "$1^{\mathrm{st}}$ Stage Mass", ...
                                "$2^{\mathrm{nd}}$ Stage Mass", ...
                                "Gross Vehicle Mass", ...
                                "Minimum Gross Vehicle Mass", ...
                                Interpreter="latex", ...
                                Location="eastoutside");


                        xlim([0 1]);
                        ylim([0 max(m_0s)/1000]);
                        hold off;
                end

                function obj = save_mass_fig(obj)
                        arguments
                                % OBJ The object itself
                                obj TwoStageLV
                        end

                        filename = sprintf("./images/s1 %s - s2 %s.jpg", obj.propellants(1).name, obj.propellants(2).name);
                        obj.mass_fig.Units = "centimeters";
                        obj.mass_fig.Position = [2 2 32 24];
                        saveas(obj.mass_fig, filename);
                end
        end
end
