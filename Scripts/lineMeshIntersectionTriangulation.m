function [isIntersecting, intersectionLength] = lineMeshIntersectionTriangulation3D(lineSegment, triangulationMatrix, vertices)
    % lineSegment: a 2x3 matrix representing the start and end points of the line segment
    % triangulationMatrix: a mx3 matrix representing the connectivity of m triangles
    % vertices: a nx3 matrix representing the coordinates of n vertices
    
    numTriangles = size(triangulationMatrix, 1);

    % Extract triangle vertices using triangulationMatrix
    triangleVertices = reshape(vertices(triangulationMatrix', :)', 3, 3, numTriangles);

    % Extract line segment start and end points
    P0 = lineSegment(1, :); % Start point of the line segment
    P1 = lineSegment(2, :); % End point of the line segment

    % Initialize result arrays
    isIntersecting = false(1, numTriangles);
    intersectionLength = zeros(1, numTriangles);

    % Compute vectors
    AB = triangleVertices(:, 2, :) - triangleVertices(:, 1, :);
    AC = triangleVertices(:, 3, :) - triangleVertices(:, 1, :);
    AP0 = P0 - triangleVertices(:, 1, :);

    for i = 1:numTriangles
        % Compute normal of the triangle
        N = cross(AB(:, :, i), AC(:, :, i), 1);

        % Check if line segment is parallel to the triangle
        if dot(N, cross(AB(:, :, i), AP0(:, :, i), 1)) ~= 0
            % Compute barycentric coordinates
            alpha = dot(N, cross(AB(:, :, i), (P1 - triangleVertices(:, 1, i)), 1)) / dot(N, N);
            beta = dot(N, cross((P0 - triangleVertices(:, 1, i)), AC(:, :, i), 1)) / dot(N, N);

            % Check if intersection is within the triangle
            if (alpha >= 0) && (beta >= 0) && (alpha + beta <= 1)
                % Calculate the intersection point
                intersectionPoint = P0 + alpha * (P1 - P0);
                % Calculate the length of the intersection segment
                intersectionLength(i) = norm(intersectionPoint - P0);
                isIntersecting(i) = true;
            end
        end
    end
end