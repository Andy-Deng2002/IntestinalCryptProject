% filepath: c:\Users\xiluo\Desktop\Crypt simulation\multi-dimensional\CheckLumpability.m
function CheckLumpability()
    % --- Configuration from  ---
    ROWS = 1;
    COLS = 12;      
    LAMBDA = 10.0;   
    ALPHA = 10.0;
    
    fprintf('--- LUMPABILITY CHECK (%dx%d, Lambda=%.1f) ---\n', ROWS, COLS, LAMBDA);
    fprintf('Definition: Sum_{y in S_{k+1}} P(x->y) must be constant for all x in S_k.\n\n');
    
    N_cells = ROWS * COLS;
    num_states = 2^N_cells;
    
    % 1. Build Exact Graph Structure
    Adj = build_hex_adjacency(ROWS, COLS);
    W = Adj / 6.0; % Probability of picking neighbor j from i (Standard Random Walk)
    
    % Build Selection Weights (Pi)
    row_w = ALPHA .^ (ROWS:-1:1);
    pi_vector = zeros(N_cells, 1);
    for r=1:ROWS, for c=1:COLS, pi_vector(sub2ind([ROWS,COLS],r,c)) = row_w(r); end, end
    pi_vector = pi_vector / sum(pi_vector);
    
    % 2. Storage: Group probabilities by mutant count k
    % P_up_stats{k+1} stores P(x -> Any State in k+1) for every x
    P_up_stats = cell(N_cells, 1); 
    
    % 3. Iterate ALL states to calculate transition sums
    for s = 0:(num_states-1)
        x = bitget(s, 1:N_cells)'; % Convert integer to vector x
        k = sum(x);
        
        if k == N_cells || k == 0, continue; end 
        
        % Calculate Denominator (Normalization constant Z_x)
        zeta = dot(x, pi_vector);
        Z = 1 + (LAMBDA - 1) * zeta;
        
        % Calculate Total Probability to go UP (Gain a mutant anywhere)
        prob_gain_total = 0;
        
        for j = 1:N_cells % Target cell (must be WT to gain)
            if x(j) == 0 
                % Sum over all neighbors i that could replace j
                nbrs = find(Adj(:, j));
                for i = nbrs'
                    if x(i) == 1 % Neighbor is Mutant
                        % P(x -> y) contribution
                        w_ij = W(i, j); 
                        p_birth = (LAMBDA * pi_vector(i)) / Z;
                        prob_gain_total = prob_gain_total + p_birth * w_ij;
                    end
                end
            end
        end
        
        % Store the sum for this specific configuration x
        P_up_stats{k+1} = [P_up_stats{k+1}; prob_gain_total];
    end
    
    % 4. specific check for k=1 to N-1
    fprintf('Mutants (k) |  Min P(up)  |  Max P(up)  |  Diff   | Lumpable?\n');
    fprintf('------------------------------------------------------------\n');
    
    is_strictly_lumpable = true;
    
    for k = 1:(N_cells-1)
        vals = P_up_stats{k+1};
        if isempty(vals), continue; end
        
        min_v = min(vals);
        max_v = max(vals);
        diff_v = max_v - min_v;
        
        % Check if difference is effectively zero (floating point tolerance)
        if diff_v > 1e-12
            status = 'NO';
            is_strictly_lumpable = false;
        else
            status = 'YES';
        end
        
        fprintf('     %2d     |   %.5f   |   %.5f   | %.1e |  %s\n', k, min_v, max_v, diff_v, status);
    end
    
    fprintf('------------------------------------------------------------\n');
    if is_strictly_lumpable
        fprintf('RESULT: YES. The system is STRICTLY LUMPABLE.\n');
        fprintf('Exact Graph P(fix) === Mean Field P(fix).\n');
    else
        fprintf('RESULT: NO. The transition rates depend on spatial configuration.\n');
        fprintf('Example: A mutant block grows slower than scattered mutants.\n');
        fprintf('However, Fixation Probability may still be identical due to Isothermal Thm.\n');
    end

end

% --- Helper: Graph Topology (Matched to Solver) ---
function Adj = build_hex_adjacency(rows, cols)
    num = rows * cols;
    u_list=[]; v_list=[]; val_list=[]; 
    for r = 1:rows
        for c = 1:cols
            u = sub2ind([rows, cols], r, c);
            % Left/Right
            v_l=sub2ind([rows,cols],r,mod(c-2,cols)+1); 
            v_r=sub2ind([rows,cols],r,mod(c,cols)+1);
            u_list=[u_list;u;u]; v_list=[v_list;v_l;v_r]; val_list=[val_list;1;1];
            
            % Vertical (if rows > 1)
            % ... (simplified for 1D case, extend if checking 2D) ...
        end
    end
    Adj = sparse(u_list, v_list, val_list, num, num);
end