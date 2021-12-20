function [coeffs, boundary_picked] = find_coefficients_var(boundary_segments, invariants, degree_poly, num_picked_clusters, cluster_size, gradient_bound, invariant_bound, solver, display)
    fprintf("\n[COEFFICIENTS]: setting up problem\n");
    t0 = tic;

    %% preparation
    syms x y % declaring the symbols
    
    monomials = monomials_gen([x;y], 0:degree_poly); % generate the monomials vector
    grad_x = diff(monomials, x);
    grad_y = diff(monomials, y);
    num_monomials = size(monomials, 1); % number of monomials
    
    boundary_clusters = cluster_neighbors(boundary_segments, cluster_size); % boundary clusters
    
    clusters_sizes = cell2mat(cellfun(@size, boundary_clusters, 'UniformOutput', false)); % size of each cluster
    clusters_picked = randsample(find(clusters_sizes(:, 1) == repmat(cluster_size, size(clusters_sizes, 1), 1)), num_picked_clusters, false); % clusters picked among those with full size
    boundary_picked = cell2mat(boundary_clusters(clusters_picked));
    num_picked = num_picked_clusters * cluster_size;
    
    num_variables =   num_monomials ... % polynomial coefficients
                    + num_picked ... % weights for picked clusters
                    + 1; % epsilon
                    % the exact definition of the variables follows this order
    num_dimension = num_variables + 1; % constant term
    
    
    %% building F
    pop.F = zeros(1, num_dimension); % zero-initialization
    
    % epsilon^2
    pop.F(1,num_dimension) = 1; pop.F(1,num_variables) = 2;
    
    for i = 1:num_picked_clusters
        % beginning and ending index of the current cluster
        idx_1 = num_monomials + (i-1) * cluster_size + 1;
        
        % sample variance var(x) = 1/n * sum{x^2} - 1/((n-1)*n) * sum{x}^2
        
        % term1 = 1/(n-1) * sum{x^2}
        term1 = zeros(cluster_size, num_dimension); term1(:, num_dimension) = 1 / (cluster_size-1);
        for j = 1:cluster_size
            term1(j, idx_1 + j - 1) = 2;
        end
        
        % sum_x = sum{x}
        sum_x = zeros(cluster_size, num_dimension); sum_x(:, num_dimension) = 1; 
        for j = 1:cluster_size
           sum_x(j, idx_1 + j - 1) = 1; 
        end
        % term2 = - 1 / n / (n-1) * sum_x ^ 2
        term2 = SBSOS_constant_multiply(- 1 / cluster_size / (cluster_size - 1), SBSOS_multiply(sum_x, sum_x));
        
        pop.F = SBSOS_add(pop.F, SBSOS_add(term1, term2));
    end
    
    %% building G
    num_ineq =   1 ... % constant term coefficient positive
               + 1 ... % epsilon > 0
               + 2 ... % coefficients have norm 1
               + 2 * num_picked ... % weight^2 - weight = 0 or weight = 1
               + 2 * num_picked ... %  vanishing on the picked samples
               + num_picked ... % non-vanishing gradient
               + 1; % invariant conditions
    
    pop.G = cell(1, num_ineq);
    
    % epsilon > 0
    pop.G{1} = zeros(1, num_dimension);
    pop.G{1}(1, num_variables) = 1; pop.G{1}(1, num_dimension) = 1;
    
    % c_1 >= 0 for identifiability
    pop.G{2} = zeros(1, num_dimension);
    pop.G{2}(1, 1) = 1; pop.G{2}(1, num_dimension) = 1;

    %     1 <= norm(coeffs) <= 1 
    % <=> norm(coeffs)^2 - 1 >= 0 & - norm(coeffs)^2 + 1 >= 0
    pop.G{3} = zeros(num_monomials+1, num_dimension);
    for i = 1:num_monomials
       pop.G{3}(i, num_dimension) = 1;
       pop.G{3}(i, i) = 2;
    end
    pop.G{3}(num_monomials+1, num_dimension) = -1;
    
    pop.G{4} = zeros(num_monomials+1, num_dimension);
    for i = 1:num_monomials
       pop.G{4}(i, num_dimension) = -1;
       pop.G{4}(i, i) = 2;
    end
    pop.G{4}(num_monomials+1, num_dimension) =  1;
    
    idx_ineq = 4;
    
    % w = 1        <=>  w - 1 >= 0   & -w + 1 >= 0,   for the middle point of each cluster
    % w^2 - w = 0  <=>  w^2 - w >= 0 & -w^2 + w >= 0, for the rest of the cluster
    idx = num_monomials;
    for i = 1:num_picked_clusters
        for j = 1:cluster_size
            idx = idx + 1;
            
            if j == ceil(cluster_size / 2)
                idx_ineq = idx_ineq + 1;
                pop.G{idx_ineq} = zeros(2, num_dimension);
                pop.G{idx_ineq}(1, idx) = 1;
                pop.G{idx_ineq}(1, num_dimension) =  1;
                pop.G{idx_ineq}(2, num_dimension) = -1;
    
                idx_ineq = idx_ineq + 1;
                pop.G{idx_ineq} = zeros(2, num_dimension);
                pop.G{idx_ineq}(1, idx) = 1;
                pop.G{idx_ineq}(1, num_dimension) = -1;
                pop.G{idx_ineq}(2, num_dimension) =  1;
            else
                idx_ineq = idx_ineq + 1;
                pop.G{idx_ineq} = zeros(2, num_dimension);
                pop.G{idx_ineq}(1, idx) = 2;
                pop.G{idx_ineq}(1, num_dimension) =  1;
                pop.G{idx_ineq}(2, idx) =  1;
                pop.G{idx_ineq}(2, num_dimension) = -1;
    
                idx_ineq = idx_ineq + 1;
                pop.G{idx_ineq} = zeros(2, num_dimension);
                pop.G{idx_ineq}(1, idx) = 2;
                pop.G{idx_ineq}(1, num_dimension) = -1;
                pop.G{idx_ineq}(2, idx) =  1;
                pop.G{idx_ineq}(2, num_dimension) =  1;
            end
        end
    end
    
    %     -e * w_i <= eval_coeff(x, y) * w_i <= e * w_i
    % <=> eval_coeff(x, y) * w_i + e * w_i >= 0 & - eval_coeff(x, y) * w_i + e * w_i >= 0
    
    % |grad|^2 > gradient_bound <=> |grad|^2 - gradient_bound > 0
    idx_w = num_monomials;
    const = zeros(1, num_dimension); const(1, num_dimension) = - gradient_bound;
    for i = 1:num_picked
        % evaluation
        idx_w = idx_w + 1;
    
        idx_ineq = idx_ineq + 1;
        pop.G{idx_ineq} = zeros(num_monomials + 1, num_dimension);
        
        evals = subs(monomials, [x, y], boundary_picked(i,:));
        for j = 1:num_monomials
            pop.G{idx_ineq}(j, j) = 1; 
            pop.G{idx_ineq}(j, idx_w) = 1;
            pop.G{idx_ineq}(j, num_dimension) = evals(j);
        end
        pop.G{idx_ineq}(num_monomials + 1, idx_w) = 1;
        pop.G{idx_ineq}(num_monomials + 1, num_dimension) = 1; 
        pop.G{idx_ineq}(num_monomials + 1, num_variables) = 1;
        
        idx_ineq = idx_ineq + 1;
        pop.G{idx_ineq} = SBSOS_negate(pop.G{idx_ineq - 1});
        pop.G{idx_ineq}(num_monomials + 1, num_dimension) = 1;
    
        % gradient
        idx_ineq = idx_ineq + 1;
        sum_x = zeros(num_monomials, num_dimension);
        sum_y = zeros(num_monomials, num_dimension);
        dx = subs(grad_x, [x, y], boundary_picked(i,:));
        dy = subs(grad_y, [x, y], boundary_picked(i,:));
        for j = 1:num_monomials
            sum_x(j, num_dimension) = dx(j); sum_x(j, j) = 1;
            sum_y(j, num_dimension) = dy(j); sum_y(j, j) = 1;
        end
        
        pop.G{idx_ineq} = SBSOS_add(SBSOS_add(SBSOS_multiply(sum_x, sum_x), SBSOS_multiply(sum_y, sum_y)), const);
    end

    %     sum c_i * c_j * inv_ij <= eps 
    % <=> eps - sum c_i * c_j * inv_ij >= 0
    idx_ineq = idx_ineq + 1;
    pop.G{idx_ineq} = zeros(num_monomials * (num_monomials + 1) / 2 + 1, num_dimension);
    k = 0;
    for i = 1:num_monomials
        k = k + 1;
        pop.G{idx_ineq}(k, num_dimension) = - invariants(i, i);
        pop.G{idx_ineq}(k, i) = 2;

        for j = (i+1):num_monomials
            k = k + 1;
            pop.G{idx_ineq}(k, num_dimension) = - 2 * invariants(i, j);
            pop.G{idx_ineq}(k, i) = 1; pop.G{idx_ineq}(k, j) = 1;
        end

        pop.G{idx_ineq}(k+1, num_dimension) = invariant_bound;
    end
    
    %% specifying parameters
    pop.n = num_variables;
    pop.k = 1;
    pop.d = 2; % 2
    pop.I = cell(1, num_picked_clusters + 1);
    pop.I{1} = [1:num_monomials, num_variables];
    for i = 1:num_picked_clusters
        pop.I{i+1} = [1:num_monomials, (num_monomials + (i-1) * cluster_size + 1):(num_monomials + i * cluster_size), num_variables];
    end
    num_I = size(pop.I, 2);

    pop.J = cell(1, num_I);
    for i = 1:num_ineq
        found = false;
        for j = 1:num_I
            indices = find(sum(abs(pop.G{i}), 1));
            if sum(ismember(indices, pop.I{j})) == size(indices, 2) - 1
                pop.J{j}(size(pop.J{j}, 2) + 1) = i;
                found = true;
                break
            end
        end
        if found == false
            disp("[COEFFICIENTS]: wrong pop.I on " + i);
            return
        end
    end
    fprintf("  => time: " + toc(t0) + "\n");

    %%
%     pop.I = {1:num_variables};
%     pop.J = {1:num_ineq};
    
    %% creating the program
    fprintf("[COEFFICIENTS]: generating SBSOS\n");
    t0 = tic;
    sdp = gendata2(pop, 'SBSOS', 'sedumi');
    fprintf("  => time: " + toc(t0) + "\n");

    %% converting and solving
    if solver == "sedumi"
        sdp.pars.eps = 0;
        sol = csol(sdp, 'sedumi');
        psol = postproc(pop, sdp, sol);
    elseif solver == "mosek"
        fprintf("[COEFFICIENTS]: converting Sedumi to Mosek\n");
        t0 = tic;
        sdp2 = convert_sedumi2mosek(sdp.At', sdp.b,sdp.C, sdp.blk);
        fprintf("  => time: " + toc(t0) + "\n");

        fprintf("[COEFFICIENTS]: solving (Mosek)\n");
        t0 = tic;
        [rcode, msg] = mosekopt('minimize echo(0)', sdp2);
        fprintf("  => time: " + toc(t0) + "\n");

        fprintf("[COEFFICIENTS]: error code: " + rcode + ', ' + msg.sol.itr.prosta + "\n");
    end

    %% recovering the results
    if solver == "sedumi"
        % outdated
        res = psol.YY{1};
        deg1 = res(2:num_dimension);
        deg2 = res((num_dimension+1):(num_dimension+num_dimension*num_variables/2));
        coeffs = res(1:num_monomials);
    elseif solver == "mosek"
        res = recovery(sdp.recy, msg.sol.itr.y);
        coeffs = res{1}(2:num_monomials+1);

        % outdated
%         deg1 = res{2}(2:num_dimension);
%         deg2 = res{2}((num_dimension+1):(num_dimension+num_dimension*num_variables/2));
    end

    %% terminate
    if display == 0, return; end
    
    %% displaying the coefficients
    disp(coeffs);
    
    %% plotting fitted curve
    scatter(boundary_picked(:, 1), boundary_picked(:, 2)); 
    % axis([-1.1 1.1 -1.1 1.1]); 
    pbaspect([1 1 1]);
    hold on; 
    % coeffs = coeffs / norm(coeffs);
    fimplicit(dot(coeffs, monomials));
    hold off;

    %% terminate
    if display == 1, return; end

    %% recovering the parameters using Cholesky factorization
    Y = zeros(num_variables);
    Y(tril(true(num_variables))) = deg2;
    Y = Y + triu(Y', 1);
    Z = Y(1:num_monomials, 1:num_monomials);
    disp(Z);
    
    %% evaluation on the sampled points
    for i = 1:num_picked
        evals = double(subs(dot(coeffs, monomials), [x, y], boundary_picked(i, :)));
        disp(evals);
    end
    
    %% gradient on the sampled points
    for i = 1:num_picked
        evals = dot(subs(grad_x, [x, y], boundary_picked(i,:)), coeffs) ^ 2 + dot(subs(grad_y, [x, y], boundary_picked(i,:)), coeffs) ^ 2;
        disp(double(evals));
    end
    
    %% weights
    % outdated
    disp(deg1((num_monomials+1):(num_monomials+num_picked)));
    
    %% epsilon
    % outdated
    disp(deg1(num_monomials+num_picked+1));
end