% cluster the boundary set into cells of the specified size
% param:
%   boundary_segments: bouncdary clusters
%   cluster_size: cluster size
% return:
%   clusters: clusters, each row contains the cell of x- and y- coordinates

function clusters = cluster_neighbors(boundary_segments, cluster_size)
    num_points = size(boundary_segments, 1);
    idx = 1; % cluster index
    
    % calculate the total number of clusters
    cnum = sum(ceil(cell2mat(cellfun(@size, boundary_segments, 'UniformOutput', false)) / cluster_size), 1);
    cnum = cnum(1);
    
    clusters = cell(cnum, 1);
    
    for i = 1:num_points % for each neighborhood cluster into csize untill finished
        num = size(boundary_segments{i}, 1);
        
        for j=1:ceil(num/cluster_size)
            clusters{idx} = boundary_segments{i} (((j-1)*cluster_size+1) : min(j*cluster_size, num), :);
            idx = idx + 1;
        end 
    end
end
