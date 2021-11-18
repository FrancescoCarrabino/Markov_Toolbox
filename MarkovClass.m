classdef MarkovClass
    % This class describes a Markov Chain from P, a non-negative squared
    % Transition Matrix, and "states", the labels for each state of the
    % Markov Chain.
    % Properties:
    %   - P: matrix. The transition matrix of a Markov Chain.
    %   - states: string array. The labels for each state of the MC.
    %   - stat_vect: double array. The stationary vector of the Markov Chain. For a
    %   non-ergodic, reducible matrix, the stat_vect is the linear
    %   combination of the stat_vect of each recurrent class. 
    %   - G: digraph. The digraph obtained from the Transition Matrix.
    %   - c_states: logical array. An array containing the classification
    %   of each state based on the connected component it belongs to and
    %   the supernode associated to that component.
    
    properties
        P;
        states;
        stat_vect;
        G;
        c_states;
    end
    
    methods
        function obj = MarkovClass(P,states)
            % Constructor Function.
            % Construct an instance of this class.
            obj.P = P;
            obj.states = states;
            obj.G = digraph(obj.P', obj.states);
            obj.stat_vect = asymptotics(obj);
            obj.c_states = classify_states(obj);
        end
        
        function ret = isreducible(obj)
            % Function that checks if the transition matrix in input is
            % reducible. It uses the theorem: 
            % Let A ∈ M^n . The following are equivalent:
            %   (a) A is irreducible.
            %   (b) (I+|A|)^(n−1) >0. 
            %   (c) (I + M(A))^(n−1) > 0.
            % From Matrix Analysis, 2nd edition, page 402.
            
            n = length(obj.P);
            I_abA = eye(n) + abs(obj.P);
            if all(I_abA^(n-1)) > 0
                ret = 0;
            else
                ret = 1;
            end
        end
        
        function ret = isergodic(obj)
            % Function that checks if matrix is ergodic, A Markov chain is 
            % ergodic if it is both irreducible and aperiodic. This 
            % condition is equivalent to the transition matrix being a 
            % primitive nonnegative matrix.
            % Ergodicity is found through the use of Wielandt's theorem.
            % From Wielandt, H. "Unzerlegbare, Nicht Negativen Matrizen." 
            % Mathematische Zeitschrift. Vol. 52, 1950, pp. 642–648.
            
            n = length(obj.states);
            m = ((n-1)^2)+1;
            
            if all(obj.P^m) > 0
                ret = 1;
            else
                ret = 0;
            end
        end
        
        function stat_vect = asymptotics(obj)
            % Power method for finding the asymptotic behaviour of the
            % Markov Chain.
            
            n = length(obj.P);
            
            if obj.isergodic()
                P_I = obj.P - eye(n);
                P_I = [P_I; ones(1, n)];
                b = zeros(1, length(P_I)-1);
                b = [b, 1]';
                stat_vect = linsolve(P_I, b);
            else
                epsilon = 0.1^21;
                v_new = ones(1,n);
                v_old = zeros(1, n, 1);
                while all(abs(v_new - v_old)) > epsilon
                    v_old = v_new;
                    v_new =  v_old*obj.P ;
                    v_new = v_new/norm(v_new);
                end
                stat_vect = v_new;
            end
         end
        
        function h = graphplot(obj)
            % Function for plotting the Markov Chain with edges coloured by
            % the corresponding transition probability and the nodes
            % highlighted in red for recurrent states and diamond shape for
            % the transition states.
            
            figure();
            colormap jet;
            h = plot(obj.G, "EdgeCData", obj.G.Edges.Weight, "LineWidth", 1.5, "MarkerSize", 10);
            cb = colorbar;
            cb.Label.String = 'Transition Probability';
            highlight(h, (obj.c_states == 1), "NodeColor", "r");
            highlight(h, (obj.c_states == 0), "Marker", "d", "NodeColor", "#4DBEEE");
            lgd = legend('Transient');
            title(lgd,'States')
            title("Markov Chain Plot")
        end
        
        function c_states = classify_states(obj)
            % Function that classifies the states of a Markov Chain in
            % "Recurrent" or "Transient", using only the stronges connected
            % components and the supernodes corresponding to the condensed
            % graph. If the supernode has outdegree == 0, then the class is
            % labeled as recurrent, else if the outdegree > 0, it is
            % labeled as transient. 
            
            c_c = conncomp(obj.G);
            c_c = F_reorder(c_c);
            cond = condensation(obj.G);
            c_c_cond = conncomp(cond);
            out_deg = centrality(cond, "outdegree");
            recurr_ids = [];
            trans_ids = [];
            c_states = ones(length(obj.P), 1);
            for i= 1:length(out_deg)
                if out_deg(i) == 0
                    recurr_ids = find(c_c == c_c_cond(i));
                else
                    
                    trans_ids = find(c_c == c_c_cond(i));
                end
            end
            for i = 1:length(recurr_ids)
                c_states(recurr_ids(i)) = 1;
            end
            for i = 1:length(trans_ids)
                c_states(trans_ids(i)) = 0;
            end
        end
        
        function h = trans_heat(obj)
            % Function which plots out the transition matrix heatmap.
            
            n = length(obj.P);
            figure;
            imagesc(obj.P);
            colormap(jet);
            colorbar;
            axis square
            h = gca;
            h.XTick = 1:n;
            h.YTick = 1:n;
            title 'Transition Matrix Heatmap';
        end
        
        function h = simulation(obj, x0, steps)
            % Function that simulates the beahviour of a user onto the
            % Markov Chain. 
            % Inputs:
            %   - obj: MarkovClass object. 
            %   - x0: double array. The initial state vector from which the
            %   user starts the simulation. Has to be inputted as row
            %   vector.
            %   - steps: integer. Number of steps for the simulation.
            
            figure();
            n = length(obj.P);
            M = zeros(steps, n);
            h = histogram(categorical(), obj.states);
            ylim([0 1]);
            title('Simulation of an Agent on the Markov Chain');
            grid
            legend("Step 0");
            if obj.isergodic()
                x_new = x0';
                for i=1:steps
                    M(i, :) = x_new';
                    h.BinCounts = x_new';
                    legend("Step " + i);
                    x_old = x_new;
                    x_new = obj.P*x_old;
                    x_new = x_new/norm(x_new);
                    pause(1);
                end
            else
                x_new = x0;
                for i=1:steps
                    M(i, :) = x_new;
                    h.BinCounts = x_new;
                    legend("Step " + i);
                    x_old = x_new;
                    x_new = x_old*obj.P;
                    x_new = x_new/norm(x_new);
                    pause(1);
                end
            end
            
            figure;
            imagesc(M);
            colormap(jet);
            colorbar;
            axis square
            hs = gca;
            hs.XTick = 1:n;
            xticklabels(obj.states);
            hs.YTick = 1:steps;
            title 'Simulation Matrix Heatmap';
        end
        
        function lc = lazy_mc(obj)
            % Function that computes a lazy Markov Chain. Takes in input a
            % Markov Chain and adds a self-loop to each state. 
            
            TM = obj.P + eye(length(obj.P));
            TM = toStochMat(TM);
            lc = MarkovClass(TM, obj.states);
        end
    end
end

