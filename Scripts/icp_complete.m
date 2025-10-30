function [aligned_nodes] = icp_complete(CP,p,iterations)
nodes_template = CP;
nodes = p;
%% Performing ICP alignment
% This is the initial alignment with no rotation.
% Two different icp approaches are used, the first includeds the faces and
% the second is just the points.
format long g
% Rotations
r.r0 = eye(3);
r.rx = rotx(90);
r.rxx = rotx(180);
r.rxxx = rotx(270);
r.ry = roty(90);
r.ryy = roty(180);
r.ryyy = roty(270);
r.rz = rotz(90);
r.rzz = rotz(180);
r.rzzz = rotz(270);
r.rxy = rotx(90) * roty(90);
r.rxyy = rotx(90) * roty(180);
r.rxyyy = rotx(90) * roty(270);
r.rxxy = rotx(180) * roty(90);
r.rxxyy = rotx(180) * roty(180);
r.rxxyyy = rotx(180) * roty(270);
r.rxxxy = rotx(270) * roty(90);
r.rxxxyy = rotx(270) * roty(180);
r.rxxxyyy = rotx(270) * roty(270);
r.rxz = rotx(90) * rotz(90);
r.rxzz = rotx(90) * rotz(180);
r.rxzzz = rotx(90) * rotz(270);
r.rxxz = rotx(180) * rotz(90);
r.rxxzz = rotx(180) * rotz(180);
r.rxxzzz = rotx(180) * rotz(270);
r.rxxxz = rotx(270) * rotz(90);
r.rxxxzz = rotx(270) * rotz(180);
r.rxxxzzz = rotx(270) * rotz(270);
r.ryx = roty(90) * rotx(90);
r.ryxx = roty(90) * rotx(180);
r.ryxxx = roty(90) * rotx(270);
r.ryyx = roty(180) * rotx(90);
r.ryyxx = roty(180) * rotx(180);
r.ryyxxx = roty(180) * rotx(270);
r.ryyyx = roty(270) * rotx(90);
r.ryyyxx = roty(270) * rotx(180);
r.ryyyxxx = roty(270) * rotx(270);
r.ryz = roty(90) * rotz(90);
r.ryzz = roty(90) * rotz(180);
r.ryzzz = roty(90) * rotz(270);
r.ryyz = roty(180) * rotz(90);
r.ryyzz = roty(180) * rotz(180);
r.ryyzzz = roty(180) * rotz(270);
r.ryyyz = roty(270) * rotz(90);
r.ryyyzz = roty(270) * rotz(180);
r.ryyyzzz = roty(270) * rotz(270);
r.rzx = rotz(90) * rotx(90);
r.rzxx = rotz(90) * rotx(180);
r.rzxxx = rotz(90) * rotx(270);
r.rzzx = rotz(180) * rotx(90);
r.rzzxx = rotz(180) * rotx(180);
r.rzzxxx = rotz(180) * rotx(270);
r.rzzzx = rotz(270) * rotx(90);
r.rzzzxx = rotz(270) * rotx(180);
r.rzzzxxx = rotz(270) * rotx(270);
r.rzy = rotz(90) * roty(90);
r.rzyy = rotz(90) * roty(180);
r.rzyyy = rotz(90) * roty(270);
r.rzzy = rotz(180) * roty(90);
r.rzzyy = rotz(180) * roty(180);
r.rzzyyy = rotz(180) * roty(270);
r.rzzzy = rotz(270) * roty(90);
r.rzzzyy = rotz(270) * roty(180);
r.rzzzyyy = rotz(270) * roty(270);
fields = fieldnames(r);

iterations_temp = 10;
for n = 1:numel(fields)
    rot = r.(fields{n}); % Access each rotation matrix using the field name
    rotnodes = nodes*rot; % Multiple nodes by rotation matrix
    [~,~,error_temp] = icp(nodes_template',rotnodes', iterations_temp,'Matching','kDtree'); % Perform small ICP
    E_short.(fields{n}) = error_temp(end); % Save the lowest error
end
% Convert the structure 'E' to a cell array for easier sorting
E_short_fields = fieldnames(E_short);
E_short_values = struct2array(E_short);
% Find the indices of the 5 smallest error values
[~, idx_smallest] = mink(E_short_values, 5);
% Get the 5 corresponding field names (rotation matrices)
smallest_fields = E_short_fields(idx_smallest);
% Rerun the loop with iterations on the 5 smallest error rotations
for i = 1:numel(smallest_fields)
    field_name = smallest_fields{i};  % Get the field name of the current rotation
    rot = r.(field_name);  % Access the corresponding rotation matrix
    rotnodes = nodes * rot;  % Multiply nodes by the rotation matrix
    [R_temp,T_temp,E_temp] = icp(nodes_template',rotnodes', iterations,'Matching','kDtree');
    [Rwr_temp,Twr_temp,Ewr_temp] = icp(nodes_template',rotnodes', iterations,'Matching','kDtree','WorstRejection',0.1);
    if E_temp(end) < Ewr_temp(end)
        R.(field_name) = R_temp;
        T.(field_name) = T_temp;
        E.(field_name) = E_temp(end);
    else
        R.(field_name) = Rwr_temp;
        T.(field_name) = Twr_temp;
        E.(field_name) = Ewr_temp(end);
    end
end

% Find the smallest error value and corresponding field name
E_values = struct2array(E);  % Convert the structure 'E' to a regular array of error values
E_fields = fieldnames(E);    % Get the list of field names from 'E'
[~, idx_smallest] = min(E_values);  % Find the index of the smallest error value
smallest_field = E_fields{idx_smallest};  % Get the corresponding field name
% Retrieve the corresponding R, T, and rotation matrix
best_R = R.(smallest_field);  % The best R matrix
best_T = T.(smallest_field);  % The best T vector
best_rotation_matrix = r.(smallest_field);  % The rotation matrix from the original structure
% Perform the final alignment calculation
aligned_nodes = (best_R * ((nodes*best_rotation_matrix)') + repmat(best_T, 1, length(nodes')))';  % Align the nodes
% % Store the results for the final transformation
% iflip = best_rotation_matrix;  % The rotation matrix used for alignment
% iR = best_R;  % The best R matrix
% iT = best_T;  % The best T vector
%% Visualize proper alignment
figure()
plot3(nodes_template(:,1),nodes_template(:,2),nodes_template(:,3),'.k')
hold on
plot3(aligned_nodes(:,1),aligned_nodes(:,2),aligned_nodes(:,3),'.b')
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal





