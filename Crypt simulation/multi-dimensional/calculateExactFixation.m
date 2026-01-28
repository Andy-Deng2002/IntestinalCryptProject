% filepath: c:\Users\xiluo\Desktop\Crypt simulation\calculate_exact_fixation.m
function [probs, u] = calculateExactFixation(rows, cols, lambda, alpha)
% CALCULATE_EXACT_FIXATION Solves the fixation probability for a multi-layer grid
% using a Birth-Death process (Selection for Proliferaion).
%
% Inputs:
%   rows   - Number of layers
%   cols   - Number of cells per ring (N)
%   lambda - Fitness advantage of mutant (>1 = beneficial)
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
    
    % --- Weight Definition (Birth/Proliferation Potential) ---
    % As per request: alpha^(N-j+1) logic mapped to rows
    % Row 1 (Bottom) gets highest weight: alpha^rows
    % Row R (Top) gets lowest weight: alpha^1
    layer_weights = alpha .^ (rows:-1:1); 
    
    % --- Generate State Space ---
    all_states = generate_states_recursive(rows, N);
    
    % --- Matrix Construction ---
    % Estimates for sparse pre-allocation
    max_nnz = num_states * (1 + 6 * rows); % Self + (Gain/Loss)*(Self+Up+Down loops)
    I_list = zeros(max_nnz, 1);
    J_list = zeros(max_nnz, 1);
    V_list = zeros(max_nnz, 1);
    b      = zeros(num_states, 1);
    k_cnt  = 0;
    
    for idx = 1:num_states
        ns = all_states(idx, :); 
        curr_idx = idx;
        
        % Boundary: Extinction (All 0)
        if all(ns == 0)
            k_cnt=k_cnt+1; I_list(k_cnt)=curr_idx; J_list(k_cnt)=curr_idx; V_list(k_cnt)=1.0;
            b(curr_idx) = 0.0;
            continue;
        end
        
        % Boundary: Fixation (All N)
        if all(ns == N)
            k_cnt=k_cnt+1; I_list(k_cnt)=curr_idx; J_list(k_cnt)=curr_idx; V_list(k_cnt)=1.0;
            b(curr_idx) = 1.0;
            continue;
        end
        
        % --- Total Weight Calculation (Denominator) ---
        % W_tot = Sum over all cells of their proliferation rates
        % Rate(WT) = 1 * weight
        % Rate(Mut) = lambda * weight
        W_tot = 0;
        for r = 1:rows
            W_tot = W_tot + layer_weights(r) * ( (N - ns(r)) + ns(r)*lambda );
        end
        
        total_rate_out = 0;
        
        % --- Transitions Loop ---
        % In Birth-Death, we iterate over the SOURCE of the birth event
        for src = 1:rows
            [deg, k_same, k_up, k_down] = get_hex_connectivity(src, rows);
            
            % =============================================================
            % EVENT TYPE 1: MUTANT EXPANSION (Source 'src' is Mutant)
            % =============================================================
            if ns(src) > 0
                % Probability that a specific Mutant in 'src' is chosen to divide
                prob_select_src_mut = ( ns(src) * layer_weights(src) * lambda ) / W_tot;
                
                % 1.A: Mutant replaces WT in SAME layer (src -> src)
                if ns(src) < N
                    target = src;
                    % Target pool size is N-1 (cannot replace self)
                    prob_replace_wt = (k_same/deg) * ((N - ns(target)) / (N-1));
                    rate = prob_select_src_mut * prob_replace_wt;
                    
                    if rate > 0
                        ns_next = ns; ns_next(target) = ns_next(target) + 1;
                        add_transition(curr_idx, get_linear_idx(ns_next, N), rate);
                    end
                end
                
                % 1.B: Mutant replaces WT in UP layer (src -> src-1)
                if src > 1
                    target = src - 1;
                    if ns(target) < N
                        % Target pool size is N (external layer)
                        prob_replace_wt = (k_up/deg) * ((N - ns(target)) / N);
                        rate = prob_select_src_mut * prob_replace_wt;
                        
                        if rate > 0
                             ns_next = ns; ns_next(target) = ns_next(target) + 1;
                             add_transition(curr_idx, get_linear_idx(ns_next, N), rate);
                        end
                    end
                end
                
                % 1.C: Mutant replaces WT in DOWN layer (src -> src+1)
                if src < rows
                    target = src + 1;
                    if ns(target) < N
                        prob_replace_wt = (k_down/deg) * ((N - ns(target)) / N);
                        rate = prob_select_src_mut * prob_replace_wt;
                        
                        if rate > 0
                             ns_next = ns; ns_next(target) = ns_next(target) + 1;
                             add_transition(curr_idx, get_linear_idx(ns_next, N), rate);
                        end
                    end
                end
            end
            
            % =============================================================
            % EVENT TYPE 2: MUTANT LOSS (Source 'src' is WT)
            % =============================================================
            if ns(src) < N
                % Probability that a specific WT in 'src' is chosen to divide
                prob_select_src_wt = ( (N - ns(src)) * layer_weights(src) ) / W_tot;
                
                % 2.A: WT replaces Mutant in SAME layer (src -> src)
                if ns(src) > 0
                    target = src;
                    prob_replace_mut = (k_same/deg) * (ns(target) / (N-1));
                    rate = prob_select_src_wt * prob_replace_mut;
                    
                    if rate > 0
                        ns_next = ns; ns_next(target) = ns_next(target) - 1;
                        add_transition(curr_idx, get_linear_idx(ns_next, N), rate);
                    end
                end
                
                % 2.B: WT replaces Mutant in UP layer (src -> src-1)
                if src > 1
                    target = src - 1;
                    if ns(target) > 0
                        prob_replace_mut = (k_up/deg) * (ns(target) / N);
                        rate = prob_select_src_wt * prob_replace_mut;
                        
                        if rate > 0
                            ns_next = ns; ns_next(target) = ns_next(target) - 1;
                            add_transition(curr_idx, get_linear_idx(ns_next, N), rate);
                        end
                    end
                end
                
                % 2.C: WT replaces Mutant in DOWN layer (src -> src+1)
                if src < rows
                    target = src + 1;
                    if ns(target) > 0
                        prob_replace_mut = (k_down/deg) * (ns(target) / N);
                        rate = prob_select_src_wt * prob_replace_mut;
                        
                        if rate > 0
                             ns_next = ns; ns_next(target) = ns_next(target) - 1;
                             add_transition(curr_idx, get_linear_idx(ns_next, N), rate);
                        end
                    end
                end
            end
            
        end % End Source Loop
        
        % Diagonal Entry: Sum of all leaving rates (Total Rate Out)
        if total_rate_out == 0
             k_cnt=k_cnt+1; I_list(k_cnt)=curr_idx; J_list(k_cnt)=curr_idx; V_list(k_cnt)=1.0;
        else
             k_cnt=k_cnt+1; I_list(k_cnt)=curr_idx; J_list(k_cnt)=curr_idx; V_list(k_cnt)=total_rate_out;
        end
    end
    
    % --- Nested Function for adding transitions to sparse list ---
    function add_transition(from_id, to_id, rate_val)
        % Add off-diagonal term: -rate
        k_cnt = k_cnt + 1;
        I_list(k_cnt) = from_id;
        J_list(k_cnt) = to_id;
        V_list(k_cnt) = -rate_val;
        
        % Add to diagonal sum (rate leaving current state)
        total_rate_out = total_rate_out + rate_val;
    end
    
    % --- Solve Linear System ---
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