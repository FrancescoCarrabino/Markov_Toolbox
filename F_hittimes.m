function tau = F_hittimes(mc, s, ret)
% Function that calculates the hitting times from node i to node s. Can
% also calculate the return time if set in the parameters. 
%   Inputs:
%   - mc: Markov Chain, required
%       Markov Chain necessary for the calculations.
%   - s: integer, required
%       The node to calculate the hitting times on. 
%   - ret: boolean, optional
%       Default: false
%       If return is true, then the hitting times are calculated from all
%       the nodes to the node s.
%       If return is false, then it is calculated the return time from node
%       s to node s.

    if nargin < 3, ret = false; end
    
    Q = mc.P;
    Q(:, s) = [];
    Q(s, :) = [];
    I_QT = eye(length(Q)) - Q';
    b = ones(length(Q), 1);
    tau_acc = linsolve(I_QT, b); % (i -> s)
        
    if ret
        if ~isreducible(mc)
            v = asymptotics(mc);
            tau = 1/v(s);
        else
            for j=1:length(tau_acc)
                tau =  mc.P(j, s) * tau_acc(j);
            end
            tau = 1 + tau;
        end        
        tau = round(tau);
        disp("Return time to node " + string(s) + ": " + string(tau));
            
    else
        tau = tau_acc;
        tau = round(tau);
        for i=1:length(tau_acc)
            if i >= s
                disp("Access time from node " +string(i+1) + " to node " + string(s) + ": " + string(tau(i)));
            else
                disp("Access time from node " +string(i) + " to node " + string(s) + ": " + string(tau(i)));
            end
        end
        
    end
end

