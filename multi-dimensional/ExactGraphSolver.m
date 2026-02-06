function [fixation_probs, full_solution] = ExactGraphSolver(rows, cols, lambda, alpha)
% EXACTGRAPHSOLVER Solves the Microscopic Spatial Moran Process accurately.

    N_cells = rows * cols;
    num_states = 2^N_cells;
    
    if N_cells > 12
        error('Grid too large for Matrix Solver.');
    end
    
    % --- 1. Setup Weights (Pi) ---
    % Standardize: scalar alpha -> vector

    row_w = alpha .^ (rows:-1:1); 

    % Construct Pi Vector (selection policy)
    % We use simple linear indexing: 1..cols is Row 1, cols+1..2*cols is Row 2.
    pi_vector = zeros(N_cells, 1);
   
    for r = 1:rows
        for c = 1:cols
            idx = sub2ind([rows, cols], r, c);
            pi_vector(idx) = row_w(r);
        end
    end
    
    % Normalize Pi so it sums to 1
    pi_vector = pi_vector / sum(pi_vector);
    
    % --- 2. Setup Connectivity  M---
    Adj = build_hex_adjacency(rows, cols);
    Degrees = sum(Adj, 2);
    M = Adj ./ Degrees; % Row stochastic
    
    % --- 3. Build Sparse Matrix A ---
    
    max_nz = num_states * (1 + N_cells); % approximate
    I_list = zeros(max_nz, 1);
    J_list = zeros(max_nz, 1);
    V_list = zeros(max_nz, 1);
    b = zeros(num_states, 1);
    ptr = 0;
    
    for s = 0:(num_states-1)
        curr_idx = s + 1;
        
        % Convert state integer to bit vector
        % x(k) = 1 if cell k is mutant
        x = bitget(s, 1:N_cells)'; 
        
        % A. Boundary Conditions
        n_mut = sum(x);
        if n_mut == 0
            % Extinction: Prob = 0
            ptr=ptr+1; I_list(ptr)=curr_idx; J_list(ptr)=curr_idx; V_list(ptr)=1.0;
            b(curr_idx) = 0.0;
            continue;
        elseif n_mut == N_cells
            % Fixation: Prob = 1
            ptr=ptr+1; I_list(ptr)=curr_idx; J_list(ptr)=curr_idx; V_list(ptr)=1.0;
            b(curr_idx) = 1.0;
            continue;
        end
        
        % B. Transient State Equation
        % rho(x) = sum_y P(x->y) rho(y) + P(x->x) rho(x)
        % rho(x) * (1 - P(x->x)) - sum_{y!=x} P(x->y) rho(y) = 0
        % (Total Rate Out) * rho(x) - sum Rate(x->y) * rho(y) = 0
        
        % Calculate Denominator Z for current state
        zeta = dot(x, pi_vector);
        Z = 1 + (lambda - 1) * zeta;
        
        total_rate_out = 0;
        
        % Iterate possible events: Cell i reproduces, replaces neighbor j
        % Only need to consider events that CHANGE state (x_i != x_j)
        
        for j = 1:N_cells
            % Neighbors of j
            nbrs = find(Adj(:, j)); 
            
            for i = nbrs'
                if x(i) ~= x(j)
                    % Rate calculation
                    m_ij = M(i, j); % Prob i places offspring in j
                    
                    if x(i) == 1
                        % Mutant Birth (i=1) -> Replaces WT (j=0) -> Gain
                        p_birth = (lambda * pi_vector(i)) / Z;
                        
                        % Calculate Target State Index
                        % Bit j flips 0->1. Value increases by 2^(j-1)
                        next_s = s + 2^(j-1);
                    else
                        % WT Birth (i=0) -> Replaces Mutant (j=1) -> Loss
                        p_birth = (pi_vector(i)) / Z;
                        
                        % Bit j flips 1->0. Value decreases by 2^(j-1)
                        next_s = s - 2^(j-1);
                    end
                    
                    rate = p_birth * m_ij;
                    
                    if rate > 0
                        % Add -rate to off-diagonal (current, next)
                        next_idx = next_s + 1;
                        ptr=ptr+1; I_list(ptr)=curr_idx; J_list(ptr)=next_idx; V_list(ptr)= -rate;
                        
                        total_rate_out = total_rate_out + rate;
                    end
                end
            end
        end
        
        % Add Diagonal (Total Rate Out)
        if total_rate_out == 0
             % Should not happen for transient states in connected graph
             ptr=ptr+1; I_list(ptr)=curr_idx; J_list(ptr)=curr_idx; V_list(ptr)=1.0;
        else
             ptr=ptr+1; I_list(ptr)=curr_idx; J_list(ptr)=curr_idx; V_list(ptr)=total_rate_out;
        end
    end
    
    % --- 4. Solve System ---
    A = sparse(I_list(1:ptr), J_list(1:ptr), V_list(1:ptr), num_states, num_states);
    full_solution = A \ b;
    
    % --- 5. Extract Result for Single Mutants ---
    fixation_probs = zeros(rows, 1);
    for r = 1:rows
        % Find index of a state with exactly 1 mutant in row r (col 1)
        idx = sub2ind([rows, cols], r, 1);
        state_int = 2^(idx-1); 
        fixation_probs(r) = full_solution(state_int + 1);
    end
end

% --- Dependency: Valid Hex Graph ---
function Adj = build_hex_adjacency(rows, cols)
    num = rows*cols;
    % Create triplets for sparse matrix accumulation
    % Upper boound of edges: 6 neighbors * num cells
    max_edges = num * 6;
    u_list = zeros(max_edges, 1);
    v_list = zeros(max_edges, 1);
    val_list = zeros(max_edges, 1);
    ptr = 0;
    
    for r = 1:rows
        for c = 1:cols
            u = sub2ind([rows, cols], r, c);
            
            % Defined Directions (Simulation Style logic)
            targets = [];
            
            % 1. Horizontal Left & Right (Always exist)
            c_left = mod(c-2, cols) + 1;
            c_right = mod(c, cols) + 1;
            
            targets = [targets; sub2ind([rows, cols], r, c_left)];
            targets = [targets; sub2ind([rows, cols], r, c_right)];
            
            % 2. Vertical Neighbors (Look UP and DOWN)
            
            % Find rows above and below
            r_up = r - 1;
            r_dn = r + 1;
            
            if mod(r, 2) == 1 % Odd Row
                % Neighbors Above (Even): c-1, c
                if r_up >= 1
                    c_ul = mod(c-2, cols) + 1; % c-1
                    c_ur = c;
                    targets = [targets; sub2ind([rows, cols], r_up, c_ul)];
                    targets = [targets; sub2ind([rows, cols], r_up, c_ur)];
                end
                
                
                % Sim: Odd Row (r) neighbors Down (r+1 Even): c, c+1
                if r_dn <= rows
                   c_dl = c;
                   c_dr = mod(c, cols) + 1; % c+1
                   targets = [targets; sub2ind([rows, cols], r_dn, c_dl)];
                   targets = [targets; sub2ind([rows, cols], r_dn, c_dr)];
                end
                
            else % Even Row
                % Neighbors Above (Odd): c, c+1 (Inverse of Odd-to-Even-Down)
                if r_up >= 1
                    c_ul = c;
                    c_ur = mod(c, cols) + 1;
                    targets = [targets; sub2ind([rows, cols], r_up, c_ul)];
                    targets = [targets; sub2ind([rows, cols], r_up, c_ur)];
                end
                
                % Even Row(r) neighbors Down(r+1 Odd): c-1, c
                if r_dn <= rows
                    c_dl = mod(c-2, cols) + 1; % c-1
                    c_dr = c;
                    targets = [targets; sub2ind([rows, cols], r_dn, c_dl)];
                    targets = [targets; sub2ind([rows, cols], r_dn, c_dr)];
                end
            end
            
            % Add to list
            for k = 1:length(targets)
                ptr = ptr + 1;
                u_list(ptr) = u;
                v_list(ptr) = targets(k);
                val_list(ptr) = 1; % Weight of direction
            end
        end
    end
    
    % Use sparse to ACCUMULATE duplicate edges (e.g., if v_left == v_right)
    Adj = sparse(u_list(1:ptr), v_list(1:ptr), val_list(1:ptr), num, num);
end