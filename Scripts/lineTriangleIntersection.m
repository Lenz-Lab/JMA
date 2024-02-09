function doesIntersect = lineTriangleIntersection(line, triangle)
    % line: 2x3 array representing the start and end points of the line
    % triangle: 3x3 array representing the vertices of the triangle

    % Calculate the normal of the triangle
    normal = cross(triangle(2,:) - triangle(1,:), triangle(3,:) - triangle(1,:));
    normal = normal / norm(normal);

    % Calculate the direction vector of the line
    direction = line(2,:) - line(1,:);
    direction = direction / norm(direction);

    % Calculate the denominator for intersection test
    denominator = dot(normal, direction);

    % If denominator is almost zero, the line and the plane are parallel
    if abs(denominator) < 1e-6
        doesIntersect = false;
        return;
    end

    % Calculate the parameter t to find the intersection point
    t = dot(triangle(1,:) - line(1,:), normal) / denominator;

    % If t is negative or greater than 1, intersection point is outside line segment
    if t < 0 || t > 1
        doesIntersect = false;
        return;
    end

    % Calculate the intersection point
    intersectionPoint = line(1,:) + t * direction;

    % Check if the intersection point is inside the triangle
    % Using barycentric coordinates
    barycentricCoords = zeros(3,1);
    for i = 1:3
        edge1 = triangle(mod(i,3)+1,:) - triangle(i,:);
        edge2 = intersectionPoint - triangle(i,:);
        barycentricCoords(i) = dot(cross(edge1, edge2), normal);
    end

    if all(barycentricCoords >= 0)
        doesIntersect = true;
    else
        doesIntersect = false;
    end
end