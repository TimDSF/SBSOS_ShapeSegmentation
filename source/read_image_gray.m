% read in image
% param:
%   file: name of file to read
% return:
%   image: image matrix, logical; 
%   boundary_segments: each cell is a connected boundary; 
%   boundary_set: boundaries as one whole matrix

function [image, boundary_set, boundary_segments] = read_image_gray(file)
    % read in rgb image => gray => binary 
    image = imbinarize(imread(file));
    scale = size(image)/2;
    
    % calculate boundaries set
    boundary_segments = bwboundaries(edge(image,'canny'),'noholes');
    
    % boundaries set normalized to zero-mean, one-bounded
    for i = 1:size(boundary_segments, 1)
       boundary_segments{i} = boundary_segments{i} - scale;
       boundary_segments{i} = boundary_segments{i} ./ scale;
    
       tmp = boundary_segments{i}(:, 1);
       boundary_segments{i}(:, 1) = boundary_segments{i}(:, 2);
       boundary_segments{i}(:, 2) = - tmp;
    end

    % boundary matrix
    boundary_set = cell2mat(boundary_segments);
end