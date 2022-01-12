function [coeffs, boundary_picked] = find_coefficients(boundary_segments, degree_poly, num_picked_clusters, cluster_size, gradient_bound, solver, display)
    fprintf("\n[COEFFICIENTS]: setting up problem\n");
    t0 = tic;

    %% sampling
    boundary_clusters = cluster_neighbors(boundary_segments, cluster_size); % boundary clusters
    
    clusters_sizes = cell2mat(cellfun(@size, boundary_clusters, 'UniformOutput', false)); % size of each cluster
    clusters_picked = randsample(find(clusters_sizes(:, 1) == repmat(cluster_size, size(clusters_sizes, 1), 1)), num_picked_clusters, false); % clusters picked among those with full size
    boundary_picked = cell2mat(boundary_clusters(clusters_picked));
    num_picked = num_picked_clusters * cluster_size;

    %% preparation
    syms x y % declaring the symbols
    
    monomials = monomials_gen([x;y], 0:degree_poly); % generate the monomials vector
    grad_x = diff(monomials, x);
    grad_y = diff(monomials, y);
    num_monomials = size(monomials, 1); % number of monomials
    
    num_variables =   num_monomials ... % polynomial coefficients
                    + num_picked ... % weights for picked clusters
                    + 1; % epsilon
                    % the exact definition of the variables follows this order
    num_dimension = num_variables + 1; % constant term
    
    %% building F
    pop.F = zeros(1, num_dimension); % zero-initialization
    
    % epsilon^2
    pop.F(1,num_dimension) = 1; pop.F(1,num_variables) = 2;
    
    %% building G
    num_ineq =   1 ... % constant term coefficient positive
               + 1 ... % epsilon > 0
               + 2 ... % coefficients have norm 1
               + 2 * num_picked ... %  vanishing on the picked samples
               + num_picked; % non-vanishing 2-norm on gradient
    
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
    
    %     -e <= eval_coeff(x, y) <= e
    % <=> eval_coeff(x, y) + e >= 0 & - eval_coeff(x, y) + e >= 0
    
    %     |grad|^2 > gradient_bound 
    % <=> |grad|^2 - gradient_bound > 0
    const = zeros(1, num_dimension); const(1, num_dimension) = - gradient_bound;
    
    for i = 1:num_picked
        % evaluation    
        idx_ineq = idx_ineq + 1;
        pop.G{idx_ineq} = zeros(num_monomials + 1, num_dimension);
        
        evals = subs(monomials, [x, y], boundary_picked(i,:));
        for j = 1:num_monomials
            pop.G{idx_ineq}(j, j) = 1; 
            pop.G{idx_ineq}(j, num_dimension) = evals(j);
        end
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
    
    %% specifying parameters
    mk = max(sum(pop.F, 2) - pop.F(:, num_dimension));
    for i = 1:num_ineq
        mk = max(mk, max(sum(pop.G{i}, 2) - pop.G{i}(:, num_dimension)));
    end
    disp("  => max_degree: " + mk);

    pop.n = num_variables;
    pop.k = ceil(mk / 2);
    pop.d = 1; % 1
    pop.I = {1:num_variables};
    pop.J = {1:num_ineq};
    
    fprintf("  => time: " + toc(t0) + "\n");

    save('pop.mat', "pop");
    
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

        fprintf("[COEFFICIENTS]: error code: " + rcode + ', ' + msg.rcodestr + "\n");
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