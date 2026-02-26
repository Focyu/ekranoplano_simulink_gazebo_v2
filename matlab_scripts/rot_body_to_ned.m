function [V_NED] = rot_body_to_ned(X)
V_body=[X(1);X(2);X(3)];
% phi=X(7); % phi
% theta=X(8); % theta
% psi=X(9); % psi

yaw=X(7); % phi
pitch=X(8); % theta
roll=X(9); % psi
% MFS
R=[cos(pitch)*cos(yaw) cos(yaw)*sin(pitch)*sin(roll)-sin(yaw)*cos(roll) sin(roll)*sin(yaw)+sin(pitch)*cos(yaw)*cos(roll);
cos(pitch)*sin(yaw) cos(yaw)*cos(roll)+sin(yaw)*sin(pitch)*sin(roll) sin(pitch)*sin(yaw)*cos(roll)-cos(yaw)*sin(roll);
-sin(pitch) sin(roll)*cos(pitch) cos(pitch)*cos(roll)];

% R=[cos(pitch)*cos(yaw) cos(yaw)*sin(pitch)*sin(roll)-sin(yaw)*cos(pitch) sin(roll)*sin(yaw)+sin(pitch)*cos(yaw)*cos(roll);
% cos(pitch)*sin(yaw) cos(yaw)*cos(roll)+sin(yaw)*sin(pitch)*sin(roll) sin(pitch)*sin(yaw)*cos(roll)-cos(yaw)*sin(roll);
% -sin(pitch) sin(roll)*cos(pitch) cos(pitch)*cos(roll)];
% 
% Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
% Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
% Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
% R=Rx'*Ry'*Rz';

% R_x = rotx(X(4))
% R_y = roty(X(5))
% R_z = rotz(X(6))
% R=R_x*R_y*R_z;

V_NED=R'*V_body;
V_NED=[V_NED(1);V_NED(2);V_NED(3)];
end

