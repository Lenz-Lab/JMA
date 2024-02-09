function intersects = lineSegmentTriangleIntersection(lineStart, lineEnd, triangleVertices)
    % lineStart: starting point of the line segment [x, y, z]
    % lineEnd: ending point of the line segment [x, y, z]
    % triangleVertices: vertices of the triangles [Nx3x3 matrix: [x1 y1 z1; x2 y2 z2; x3 y3 z3]]
    
    %%
    % Vector representing the line segment
    lineVec = lineEnd - lineStart;

    intersects = false;%(size(triangleVertices, 1), 1); % Initialize intersection flags
    
    for k = 1:size(triangleVertices, 1) % Iterate through each triangle
        triangle = triangleVertices(k, :, :);
        for i = 1:3 % Iterate through each triangle edge
            % Triangle edge vertices
            v0 = triangle(:,:,i);
            v1 = triangle(:,:,mod(i,3)+1); % Next vertex
           
            % Calculate vectors for edges of the triangle
            edge1 = v1 - v0;
            edge2 = triangle(:,:,mod(i+1,3)+1) - v0;
            
            % Calculate normal to the triangle
            normal = cross(edge1, edge2);
            
            % Check if the line and triangle are parallel
            if dot(normal, lineVec) == 0
                continue; % No intersection
            end
            
            % Calculate distance from the line start to the plane of the triangle
            t = dot(normal, v0 - lineStart) / dot(normal, lineVec);
            
            if t < 0 || t > 1
                continue; % Intersection is outside the line segment
            end
            
            % Calculate intersection point
            intersectionPoint = lineStart + t * lineVec;
            
            % Check if the intersection point lies within the triangle
            u = dot(intersectionPoint - v0, edge1) / dot(edge1, edge1);
            v = dot(intersectionPoint - v0, edge2) / dot(edge2, edge2);
            
            if u >= 0 && v >= 0 && u + v <= 1
                intersects          = true; % Intersection found
                intersectionPoint   = lineStart + t * lineVec;
                intersectionDist    = t;
                break; % No need to check other edges
            end
        end
    end
end