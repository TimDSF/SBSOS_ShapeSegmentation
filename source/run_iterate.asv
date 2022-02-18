%% initialization
addpaths;
clear; clc;
rng(1);

%% parameters
% input_file = 'images/aircraft_1.png';
% input_file = 'images/lip.png';
input_file = 'images/circle.png';
degree_poly = 4; % 6
num_picked_clusters = 2;
cluster_size = 10;
cluster_sparsity = 2;
gradient_bound = 0.15; % 0.05
invariant_bound = 0.50; % 1
hausdorff_lowerbound = 0.1;
hausdorff_upperbound = 0.5;
iterations = 4;
solver = "mosek";

enable_invariant = true;

dists = zeros(1, iterations);

%% training data
[image, boundary_set, boundary_segments] = read_image_gray(input_file); % read image
syms x y;
monomials = monomials_gen([x;y], 0:degree_poly); % generate the monomials vector
num_monomials = size(monomials, 1);

%% 
if isfile("_samples.mat")
    load("_samples.mat", "samples")
else
    samples = cell(1, 2);
    save("_samples.mat", "samples");
end

%% curve-fitting

for i = 1:iterations
    fprintf("\n--------------------------------");
    fprintf("\n[MAIN]: running on example " + i + "\n");

    if size(samples{2}, 2) >= 3 && enable_invariant
        invariants = find_invariants(samples, num_monomials, invariant_bound, 0);
        [coeffs, boundary_picked] = find_coefficients_var(boundary_segments, invariants, degree_poly, num_picked_clusters, cluster_size, cluster_sparsity, gradient_bound, invariant_bound, solver, 0);
        fprintf("\n[MAIN]: Invariant: " + trace(coeffs * coeffs' * invariants) + "\n");
    else
        [coeffs, boundary_picked] = find_coefficients(boundary_segments, degree_poly, num_picked_clusters, cluster_size, cluster_sparsity, gradient_bound,solver, 0);
    end

    fprintf("  => norm: " + norm(coeffs)+"\n");
    coeffs = coeffs / norm(coeffs);
    
    %%
    roots = function_roots(dot(coeffs, monomials));
    if size(roots, 1) == 0
        dists(i) = 100;
    else
        dists(i) = HausdorffDist(boundary_set, roots);
    end
    fprintf("\n[MAIN]: Hausdoff Distance: " + dists(i));

    %% plot results
    hold on;
    scatter(boundary_set(:, 1), boundary_set(:, 2), '.');
    scatter(boundary_picked(:, 1), boundary_picked(:, 2), 'o');
    fimplicit(dot(coeffs, monomials));
%     scatter(roots(:, 1), roots(:, 2), '.');
    axis equal;
    pbaspect([1 1 1]);
    hold off;

    %%
    if dists(i) < hausdorff_lowerbound
        fprintf("\n[MAIN]: positive example\n");
        samples{1}(:, size(samples{1}, 2)+1) = coeffs;
        savefig("./res/" + i + "+");
    elseif dists(i) == 100
        fprintf("\n[MAIN]: void example\n");
        samples{2}(:, size(samples{2}, 2)+1) = coeffs;
        savefig("./res/" + i + "x");
    elseif dists(i) > hausdorff_upperbound
        fprintf("\n[MAIN]: negative example\n");
        samples{2}(:, size(samples{2}, 2)+1) = coeffs;
        savefig("./res/" + i + "-");
    else
        fprintf("\n[MAIN]: undetermined example\n");
        savefig("./res/" + i);
    end
    fprintf("\n[MAIN]: samples count: [%d / %d]\n", size(samples{1}, 2), size(samples{2}, 2));

    save("_samples.mat", "samples");
end

%% save results
close;

hold on;
plot(dists);
xlim([0 iterations]); ylim([0 5]); pbaspect([2 1 1]);
plot(hausdorff_lowerbound * ones(1, iterations));
plot(hausdorff_upperbound * ones(1, iterations));
hold off;
savefig("res/dists");

return;

%% helper: 1
for i = 1:size(samples{1}, 2)
    hold on;
    f = fimplicit(dot(samples{1}(:, i), monomials)); f.Color = 'b';
    f = scatter(boundary_set(:, 1), boundary_set(:, 2), '.', 'black');
    pbaspect([1 1 1]);
    axis([-1.1 1.1 -1.1 1.1]); 
    hold off;
    pause;
    close;
end

for i = 1:size(samples{2}, 2)
    hold on;
    f = fimplicit(dot(samples{2}(:, i), monomials)); f.Color = 'r';
    scatter(boundary_set(:, 1), boundary_set(:, 2), '.', 'black');
    pbaspect([1 1 1]);
    axis([-1.1 1.1 -1.1 1.1]); 
    hold off;
    pause;
    close;
end

%% helper: 2
load("_samples.mat", "samples");
samples{1} = samples{1}(:, 1:3);
samples{2} = ones(num_monomials, 0);
save("_samples.mat", "samples");
