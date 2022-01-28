function [X] = find_invariants(samples, num_monomials, eps, display)    
    %% variables
    num_positives = size(samples{1}, 2);
    num_negatives = size(samples{2}, 2);
    num_samples = num_positives + num_negatives;
    samples = [samples{1}, samples{2}];
    
    %% Scalars
    prob.c = [];
    prob.a = sparse([], [], [], num_samples, 0);
    prob.bardim = [num_monomials];
    
    %% Objective
    prob.barc.subj = [];
    prob.barc.subk = [];
    prob.barc.subl = []; 
    prob.barc.val = [];
    
    %% Constraints
    prob.blc = [zeros(1, num_positives), ones(1, num_negatives) / eps];
    prob.buc = [ ones(1, num_positives) * eps,  inf(1, num_negatives)];
    
    prob.bara.subi = [];
    for i = 1:num_samples
        A = samples(:, i) * samples(:, i)';
        [r, c, v] = find(tril(A));
    
        idx1 = size(prob.bara.subi, 2); idx2 = idx1 + size(r, 1); idx1 = idx1 + 1;
    
        prob.bara.subi(idx1:idx2) = i;
        prob.bara.subk(idx1:idx2) = r;
        prob.bara.subl(idx1:idx2) = c;
        prob.bara.val (idx1:idx2) = v;
    end
    
    prob.bara.subj = ones(size(prob.bara.subi));
    
    %% Solve with log output
    fprintf("\n[INVARIANTS]: solving (Mosek)\n");
    t0 = tic;
    [rcode, res] = mosekopt('minimize echo(0)', prob);
    fprintf("  => time: " + toc(t0) + "\n");

    fprintf("[INVARIANTS]: error code: " + rcode + ', ' + res.sol.itr.prosta + "\n");
    
    %% Retrieve the result assuming primal and dual feasible
    X = zeros(num_monomials);
    X(tril(true(num_monomials))) = res.sol.itr.barx(1 : ((num_monomials+1) * num_monomials / 2));
    X = X + triu(X',1);
    
    %% display solutions
    if display == 0, return; end

    disp(res.sol.itr.prosta);
    disp(X);
    for i = 1:num_samples
        disp(trace(X * (samples(:, i) * samples(:, i)')));
    %     fimplicit(dot(samples(:, i), monomials));
    %     pause;
    end
end