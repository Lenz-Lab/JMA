function [isIntersecting, intersectionLength] = lineMeshIntersectionOptimized(lineSegment, triangleVertices)
    % lineSegment: a 2x3 matrix representing the start and end points of the line segment
    % triangleVertices: a nx9 matrix representing the vertices of n triangles
    
    numTriangles = size(triangleVertices, 1);

    % Reshape vertices matrix for consistency
    triangleVertices = reshape(triangleVertices', 3, 3, numTriangles);

    % Expand line segment to match the size of the triangleVertices
    P0 = repmat(lineSegment(:, 1), 1, numTriangles);
    P1 = repmat(lineSegment(:, 2), 1, numTriangles);

    % Compute vectors
    AB = triangleVertices(:, 2, :) - triangleVertices(:, 1, :);
    AC = triangleVertices(:, 3, :) - triangleVertices(:, 1, :);
    AP0 = P0 - triangleVertices(:, 1, :);
    AP1 = P1 - triangleVertices(:, 1, :);

    % Compute normal of the triangle
    N = cross(AB, AC);

    % Check if line segment is parallel to the triangles
    parallelMask = (dot(N, cross(AB, AP0)) == 0);

    % Initialize result arrays
    isIntersecting = false(1, numTriangles);
    intersectionLength = zeros(1, numTriangles);

    % Check for intersection only for non-parallel triangles
    nonParallelIndices = find(~parallelMask);
    for i = nonParallelIndices
        % Compute barycentric coordinates
        alpha = dot(N(:, i), cross(AB(:, i), AP1(:, i))) / dot(N(:, i), N(:, i));
        beta = dot(N(:, i), cross(AP0(:, i), AC(:, i))) / dot(N(:, i), N(:, i));

        % Check if intersection is within the triangle
        isIntersecting(i) = (alpha >= 0) && (beta >= 0) && (alpha + beta <= 1);

        % If there's an intersection, calculate the length
        if isIntersecting(i)
            % Calculate the intersection point
            intersectionPoint = P0(:, i) + alpha * (P1(:, i) - P0(:, i));
            % Calculate the length of the intersection segment
            intersectionLength(i) = norm(intersectionPoint - P0(:, i));
        end
    end
end