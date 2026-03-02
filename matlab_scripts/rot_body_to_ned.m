
function [V_NED] = rot_body_to_ned(X)
    V_body = [X(1); X(2); X(3)];
    
    roll = X(7);  % phi (corregido)
    pitch = X(8); % theta
    yaw = X(9);   % psi (corregido)
    
    % Matriz de Rotación de BODY a NED
    R = [cos(pitch)*cos(yaw), sin(roll)*sin(pitch)*cos(yaw)-cos(roll)*sin(yaw), cos(roll)*sin(pitch)*cos(yaw)+sin(roll)*sin(yaw);
         cos(pitch)*sin(yaw), sin(roll)*sin(pitch)*sin(yaw)+cos(roll)*cos(yaw), cos(roll)*sin(pitch)*sin(yaw)-sin(roll)*cos(yaw);
         -sin(pitch),         sin(roll)*cos(pitch),                             cos(roll)*cos(pitch)];
    
    V_NED = R * V_body; % Multiplicación directa, SIN transponer
end
