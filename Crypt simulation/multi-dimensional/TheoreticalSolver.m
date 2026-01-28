% filepath: c:\Users\xiluo\Desktop\Crypt simulation\calculate_exact_fixation.m
function [probs, u] = TheoreticalSolver(rows, cols, lambda, alpha)
% CALCULATE_EXACT_FIXATION Solves the fixation probability for a multi-layer grid.
%
% Inputs:
%   rows   - Number of layers (e.g., 3)
%   cols   - Number of cells per ring (N, e.g., 8 or 16)
%   lambda - Relative fitness cost (<1 is beneficial)
%   alpha  - Fitness stratification factor (e.g., 2.0)
%
% Outputs:
%   probs  - [rows x 1] vector containing Fixation Prob starting with 1 mutant in each layer.
%   u      - Full solution vector for all states (linear index).

    N = cols;
    dim = N + 1;
    num_states = dim^rows;
    
    if num_states > 50000 
        warning('State space is very large (%d). Calculation may use high memory.', num_states);
    end
    
    % Weight(i) = alpha^i
    layer_weights = alpha .^ (rows:-1:1); 
    
    % --- Generate State Space ---
    all_states = generate_states_recursive(rows, N);
    
    % --- Matrix Construction ---
    max_nnz = num_states * (1 + 2 * rows); 
    I_list = zeros(max_nnz, 1);
    J_list = zeros(max_nnz, 1);
    V_list = zeros(max_nnz, 1);
    b      = zeros(num_states, 1);
    k_cnt  = 0;
    
    for idx = 1:num_states
        ns = all_states(idx, :); 
        curr_idx = idx;
        
        % Boundary: Extinction
        if all(ns == 0)
            k_cnt=k_cnt+1; I_list(k_cnt)=curr_idx; J_list(k_cnt)=curr_idx; V_list(k_cnt)=1.0;
            b(curr_idx) = 0.0;
            continue;
        end
        
        % Boundary: Fixation
        if all(ns == N)
            k_cnt=k_cnt+1; I_list(k_cnt)=curr_idx; J_list(k_cnt)=curr_idx; V_list(k_cnt)=1.0;
            b(curr_idx) = 1.0;
            continue;
        end
        
        % --- Transition Rates --- (just follow the derived expressions)
        W_tot = 0;
        for r = 1:rows
            W_tot = W_tot + layer_weights(r) * ( (N - ns(r)) + ns(r)*(lambda) );
        end
        
        total_rate_out = 0;
        
        for r = 1:rows
            [deg, k_same, k_up, k_down] = get_hex_connectivity(r, rows);
            
            % Gain (n_r -> n_r + 1)
            if ns(r) < N
                rate_proliferate = ( ns(r) * layer_weights(r) * lambda ) / W_tot;
                
                prob_replace = 0;
                if ns(r) > 0,      prob_replace = prob_replace + (k_same/deg) * ((N-ns(r))/(N-1)); end
                if r > 1,          prob_replace = prob_replace + (k_up/deg) * ((N-ns(r-1))/N); end
                if r < rows,       prob_replace = prob_replace + (k_down/deg) * ((N-ns(r+1))/N); end
                
                rate_gain = rate_proliferate * prob_replace;
                
                if rate_gain > 0
                    ns_next = ns; ns_next(r) = ns_next(r) + 1;
                    next_id = get_linear_idx(ns_next, N);
                    
                    k_cnt=k_cnt+1; 
                    I_list(k_cnt)=curr_idx; J_list(k_cnt)=next_id; V_list(k_cnt)= -rate_gain;
                    total_rate_out = total_rate_out + rate_gain;
                end
            end
            
            % Loss (n_r -> n_r - 1)
            if ns(r) > 0
                rate_proliferate = ( (N - ns(r)) * layer_weights(r) ) / W_tot;
                
                prob_replace_wt = 0;
                if ns(r) < N,      prob_replace_wt = prob_replace_wt + (k_same/deg) * (ns(r)/(N-1)); end
                if r > 1,          prob_replace_wt = prob_replace_wt + (k_up/deg) * (ns(r-1)/N); end
                if r < rows,       prob_replace_wt = prob_replace_wt + (k_down/deg) * (ns(r+1)/N); end
                
                rate_loss = rate_proliferate * prob_replace_wt;
                
                if rate_loss > 0
                    ns_next = ns; ns_next(r) = ns_next(r) - 1;
                    next_id = get_linear_idx(ns_next, N);
                    
                    k_cnt=k_cnt+1; 
                    I_list(k_cnt)=curr_idx; J_list(k_cnt)=next_id; V_list(k_cnt)= -rate_loss;
                    total_rate_out = total_rate_out + rate_loss;
                end
            end
        end 
        
        % Diagonal
        if total_rate_out == 0
             k_cnt=k_cnt+1; I_list(k_cnt)=curr_idx; J_list(k_cnt)=curr_idx; V_list(k_cnt)=1.0;
        else
             k_cnt=k_cnt+1; I_list(k_cnt)=curr_idx; J_list(k_cnt)=curr_idx; V_list(k_cnt)=total_rate_out;
        end
    end
    
    % Solve
    I_list = I_list(1:k_cnt);
    J_list = J_list(1:k_cnt);
    V_list = V_list(1:k_cnt);
    A = sparse(I_list, J_list, V_list, num_states, num_states);
    u = A \ b;
    
    % Extract result for single mutants
    probs = zeros(rows, 1);
    for r = 1:rows
        s_vec = zeros(1, rows);
        s_vec(r) = 1;
        idx = get_linear_idx(s_vec, N);
        probs(r) = u(idx);
    end
end

% --- Helpers ---
function states = generate_states_recursive(num_rows, N)
    if num_rows == 1
        states = (0:N)';
    else
        sub_states = generate_states_recursive(num_rows - 1, N);
        dim_sub = size(sub_states, 1);
        dim_total = dim_sub * (N+1);
        states = zeros(dim_total, num_rows);
        row_idx = 0;
        for val = 0:N
            states(row_idx+1 : row_idx+dim_sub, 1) = val;
            states(row_idx+1 : row_idx+dim_sub, 2:end) = sub_states;
            row_idx = row_idx + dim_sub;
        end
    end
end

function idx = get_linear_idx(vec, N)
    dim = N + 1;
    idx = 1;
    R = length(vec);
    multiplier = 1;
    for i = R:-1:1
        idx = idx + vec(i) * multiplier;
        multiplier = multiplier * dim;
    end
end

function [deg, k_same, k_up, k_down] = get_hex_connectivity(r, total_rows)
    k_same = 2; 
    k_up = 0; if r > 1, k_up = 2; end
    k_down = 0; if r < total_rows, k_down = 2; end
    deg = k_same + k_up + k_down;
end