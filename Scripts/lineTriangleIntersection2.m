function intersection = lineTriangleIntersection2(p1, p2, v1, v2, v3)
    % Define the edges of the triangle
    e1 = v2 - v1;
    e2 = v3 - v1;
    
    % Compute normal vector of the triangle
    N = cross(e1, e2);
    
    % Compute determinant
    det = dot(p2 - p1, N);
    
    % Check if the line is parallel to the triangle
    if abs(det) < eps
        intersection = false;
        return;
    end
    
    % Compute direction vector of the line
    dir = p2 - p1;
    
    % Compute parameter t
    t = dot(v1 - p1, N) / det;
    
    % Check if the intersection point is outside the line segment
    if t < 0 || t > 1
        intersection = false;
        return;
    end
    
    % Compute intersection point
    intersection_point = p1 + t * dir;
    
    % Check if intersection point is inside the triangle
    edge1 = cross(e1, intersection_point - v1);
    edge2 = cross(e2, intersection_point - v2);
    edge3 = cross(v3 - v1, intersection_point - v3);
    
    if dot(N, edge1) >= 0 && dot(N, edge2) >= 0 && dot(N, edge3) >= 0
        intersection = true;
    else
        intersection = false;
    end
end