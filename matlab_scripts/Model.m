function [XDOT] = Model(X,U)

%-----------STATE AND CONTROL VECTORS-------------
x1 = X(1); %u
x2 = X(2); %v
x3 = X(3); %w
x4 = X(4); %p
x5 = X(5); %q
x6 = X(6); %r
x7 = X(7); %phi
x8 = X(8); %theta
x9 = X(9); %psi
x10 = X(10); % x_NED
x11 = X(11); % y_NED
x12 = X(12); % z_NED

u1 = U(1); %d_A (Aileron)
u2 = U(2); %d_E (Elevator)
u3 = U(3); %d_R (Rudder)
u4 = U(4); %d_th1 (throttle 1)
u5 = U(5); %d_th2 (throttle 2)

%-----------SYSTEM PARAMETERS------------------------------
m = 112; % Aircraft mass (Kg)
vc = 28; % cruise speed
bw = 5.02; % Span (Wing)
bh = 2.74; % Span (HTP)
cbar = 0.646; % Mean Aerodynamic Chord (m)
lt = 24.8;  % Distance of AC of tail and body (m)
Sw = 3.334;    % Wing planform area (m^2)
Sh = 1.128;    % Tail planform area (m^2)
ARw = 7.5; % Apect Ratio (Wing) 
ARh = 6.7; % Apect Ratio (HTP) 
TRw = 0.4; % Taper Ratio (Wing) 
TRh = 0.5; % Taper Ratio (HTP) 

Xcg = 0;  % x position of CG in Fm (m)
Ycg = 0;  % y position of CG in Fm (m)
Zcg = 0;  % z position of CG in Fm (m)

Xac = 0.097;  % x position of Aerodynamic Center in Fm (m)
Yac = 0;      % y position of Aerodynamic Center in Fm (m)
Zac = -0.06;  % z position of Aerodynamic Center in Fm (m)

%Engine constants

Xapt1 = -0.519; % x position of engine 1 force in Fm (m)
Yapt1 = 0.6;    % y position of engine 1 force in Fm (m)
Zapt1 = 0.285;  % z position of engine 1 force in Fm (m)

Xapt2 = -0.519;  % x position of engine 2 force in Fm (m)
Yapt2 = -0.6;    % y position of engine 2 force in Fm (m)
Zapt2 = 0.285;   % z position of engine 2 force in Fm (m)

%Other constants
rho = 1.225;
g = 9.81;
e0w = 0.9; % Oswald's efficiency Factor (Wing)
e0h = 0.9; % Oswald's efficiency Factor (HTP)
alpha_0w = -3.75*pi/180; % Zero Lift Angle of Attack (Wing)
alpha_0h = -4.25*pi/180; % Zero Lift Angle of Attack (HTP)
iw = 2.0 * pi/180; % Angle of incidence (Wing)
ih = 0.5 * pi/180; % Angle of Incidence (HTP)
eps = 0; % Downwash Angle
zw = 0.363; % Wing Offset
zh = 0.72; % HTP Offset
CD_0 = 0.0306; % Zero Lift Drag Coefficient
Tp = 0.125*10^(-3); % motor PWM switching frequency
%-----------CONTROL LIMITS SATURATION------------------------------

% aileron
u1min = -20*pi/180;
u1max =  20*pi/180;
% elevator
u2min = -20*pi/180;
u2max =  20*pi/180;

% rudder
u3min = -15*pi/180;
u3max =  15*pi/180;

% engines 
u4min =  0;
u4max =  Tp;

u5min =  0;
u5max = Tp;
if(u1>u1max)
    u1=u1max;
elseif(u1<u1min)
    u1=u1min;
end

if(u2>u2max)
    u2=u2max;
elseif(u2<u2min)
    u2=u2min;
end

if(u3>u3max)
    u3=u3max;
elseif(u3<u3min)
    u3=u3min;
end

if(u4>u4max)
    u4=u4max;
elseif(u4<u4min)
    u4=u4min;
end

if(u5>u5max)
    u5=u5max;
elseif(u5<u5min)
    u5=u5min;
end

%-------------INTERMEDIATE VARIABLES-----------------------------

% Va = sqrt(x1^2+ x2^2 +x3^2);   %true airspeed
% 
% alpha = atan2(x3,x1); % alpha
% beta = asin(x2/Va);   %beta
%% When Va = Vg + Vw
u = x1; %
v = x2; 
w = x3;
phi = x7; 
theta = x8; 
psi = x9;
c_phi = cos(phi);
s_phi = sin(phi);
c_theta = cos(theta);
s_theta = sin(theta);
c_psi = cos(psi);
s_psi = sin(psi);
% rotation matrix
Rb_v = [c_theta*c_psi c_theta*s_psi -s_theta;
        s_phi*s_theta*c_psi-c_phi*s_psi s_phi*s_theta*s_psi+c_phi*c_psi s_phi*c_theta;
        c_phi*s_theta*c_psi+s_phi*s_psi c_phi*s_theta*s_psi-s_phi*c_psi c_phi*c_theta];

Vb_g = [u;v;w];
% steady-state component of wind
w_ns = 0;
w_es = 0;
w_ds = 0;
% gust component of wind
u_wg = 0;
v_wg = 0;
w_wg = 0;
Vb_w = Rb_v*[w_ns;w_es;w_ds]+[u_wg;v_wg;w_wg];
% Va = Vg-Vw; 
% Va_b = [ur;vr;wr] = [u-u_w;v-v_w;w-w_w];
Va_b = [Vb_g(1)-Vb_w(1);Vb_g(2)-Vb_w(2);Vb_g(3)-Vb_w(3)];

u_r = Va_b(1);
v_r = Va_b(2);
w_r = Va_b(3);

Va = sqrt(u_r^2+v_r^2+w_r^2);
alpha = atan2(w_r,u_r);
beta = asin(v_r/Va);
%%
h = -x12;
hw = h-zw; % wing below CG -> h-offset
hh = h+zh; % tail above CG -> h+offset

Q = 0.5*rho*Va^2;     %dynamic pressure

wbe_b = [x4;x5;x6];
V_b = [x1;x2;x3];

%--------------Aerodynamic Force Coefficients---------------------

% Lift
mu_Lw = 1 + (1 - 2.25 * (TRw^0.00273 - 0.997)*((ARw^0.717)+13.6))*(288 * abs(hw/bw)^0.787 * exp(-9.14*(abs(hw/bw)^0.327)))/(ARw^0.882);
mu_Lh = 1 + (1 - 2.25 * (TRh^0.00273 - 0.997)*((ARh^0.717)+13.6))*(288 * abs(hh/bh)^0.787 * exp(-9.14*(abs(hh/bh)^0.327)))/(ARh^0.882);

C_L_e = -0.00141*u2^2+0.0307*u2;

CL_w_OGE = 2*pi*(ARw/(ARw+2))*(alpha-alpha_0w+iw);
CL_h_OGE = 2*pi*(ARh/(ARh+2))*(alpha-alpha_0h+ih+eps);

CL_w_IGE = CL_w_OGE*mu_Lw;
CL_h_IGE = (CL_h_OGE + C_L_e)*mu_Lh;

% Lh = 0.5*rho*vc^2*CL_h_IGE*Sh;
% Lw = 0.5*rho*vc^2*CL_w_IGE*Sw;
Lh = 0.5*rho*Va^2*CL_h_IGE*Sh;
Lw = 0.5*rho*Va^2*CL_w_IGE*Sw;
Ltot = Lw + Lh;

% Drag
mu_Dw = 1 - (1 - 0.157*max(0, (TRw^0.757 - 0.373))*max(0,(ARw^0.417 - 1.27)))*exp(-4.74*max(0,abs(hw/bw)^0.814))-abs(hw/bw)^2*exp(-3.88*max(0,abs(hw/bw)^0.758));
mu_Dh = 1 - (1 - 0.157*max(0, (TRh^0.757 - 0.373))*max(0,(ARh^0.417 - 1.27)))*exp(-4.74*max(0,abs(hh/bh)^0.814))-abs(hh/bh)^2*exp(-3.88*max(0,abs(hh/bh)^0.758));

C_D_e = -0.0000108*u2^2+0.000715*u2;

CD_iw_IGE = CL_w_IGE^2*mu_Dw/(pi*e0w*ARw);
CD_ih_IGE = CL_h_IGE^2*mu_Dh/(pi*e0h*ARh);

% Dtot = 0.5*rho*vc^2*(CD_0*Sw + CD_iw_IGE*Sw + CD_ih_IGE*Sh);
Dtot = 0.5*rho*Va^2*(CD_0*Sw + CD_iw_IGE*Sw + CD_ih_IGE*Sh + C_D_e*Sh);

% Side forces
CQ = -0.019*beta*180/pi; % steady-state

LD_ratio = Ltot/Dtot;
%%
%CL_wingbody
% if alpha<=alpha_switch
%     CL_wb = n*(alpha -alpha_L0);
% else
%     CL_wb = a3*alpha^3 + a2*alpha^2 + a1*alpha + a0;
% end
% 
% %CL_tail
% epsilon = depsda*(alpha - alpha_L0);
% alpha_t = alpha - epsilon + u2 + 1.3*x5*lt/Va;
% CL_t = 3.1*(St/Sw)*alpha_t;
% 
% %Total Lift coefficient
% CL = CL_wb + CL_t;
% 
% %Total Drag force
% CD = 0.13 + 0.07*(5.5*alpha + 0.654)^2;
% 
% %Side forces
% CY = -1.6*beta + 0.24*u3;

%------------------------4. Dimensional Aero Forces--------------------

% Actual dimensional forces in F_s (Stability axis)

% FA_s = [-CD*Q*Sw;
%          CQ*Q*Sw;
%         -CL*Q*Sw];
FA_s = [-Dtot;
         CQ*Q*Sw;
        -Ltot];

% Rotate forces to F_b (body axis)

C_bs = [cos(alpha) 0 -sin(alpha);0 1 0;sin(alpha) 0 cos(alpha)];   % rotation matrix

% Sequential Rotations
% Rotation about the z-axis (sideslip angle, beta): 
% R_z = [cos(beta), sin(beta), 0; 
%        -sin(beta), cos(beta), 0; 
%        0, 0, 1];
% % Rotation about the y-axis (angle of attack, alpha):
% R_y = [cos(alpha), 0, -sin(alpha); 
%        0, 1, 0;  
%        sin(alpha), 0, cos(alpha)];
% Raero→body = Ry(α)⋅Rz(β)
% C_bs = R_y*R_z;
FA_b = C_bs*FA_s;

%------------------------5. Aero Moment Coefficient about AC--------------
%calculate moments in F_b

% eta11 = -1.4*beta;
% eta21 = -0.59 - (3.1*(St*lt)/(Sw*cbar))*(alpha - epsilon);
% eta31 = (1-alpha*(180/(15*pi)))*beta;
% 
% eta = [eta11;eta21;eta31];
% 
% dCMdx = (cbar/Va)*[-11 0 5;0 (-4.03*(St*lt^2)/(Sw*cbar^2)) 0;1.7 0 -11.5];
% 
% dCMdu = [-0.6 0 0.22;0 (-3.1*(St*lt)/(Sw*cbar)) 0;0 0 -0.63];
% 
% % CM = [Cl;Cm;Cn] about Aero center in F_b
% CMac_b = eta + dCMdx*wbe_b + dCMdu*[u1;u2;u3];

% % steady-state version Sim
% Cl = -0.0005*beta*180/pi;
% Cm = -0.02*alpha*180/pi;
% Cn = -0.002*beta*180/pi;
%------------------ Cálculo de Torques Aerodinámicos ------------------
% Parámetros de estabilidad longitudinal (usando Radianes para la dinámica pura)
Cm_alpha = -1.14; % Equivalente aproximado a -0.02 * 180/pi
Cm_q = -5.0;      % Amortiguación del cabeceo 
Cm_de = -3.0;     % Autoridad del elevador (mantenemos esto en radianes también para consistencia)

% Factor de Momento por Efecto Suelo (Pitch-down moment debido al AC shift)
% A medida que hw/bw se acerca a 0, el momento de cabeceo negativo aumenta.
Cm_h = 0.05;      % Constante empírica del desplazamiento del AC
Delta_Cm_IGE = -Cm_h * exp(-4.0 * abs(hw/bw)); % Momento inducido por el suelo

Cl = -0.0005*beta*180/pi; % Eje de alabeo (Roll)

% Nuevo Cm: Momento por alfa + Momento por IGE + Momento de velocidad de cabeceo (q) + Elevador
Cm = Cm_alpha*alpha + Delta_Cm_IGE + Cm_q*(x5*cbar/(2*Va)) + Cm_de*u2; 

Cn = -0.002*beta*180/pi; % Eje de guiñada (Yaw)

CMac_b = [Cl;Cm;Cn];% steady-state
%-------------------------6. Aero moment about CG--------------------------
MAcg_b = CMac_b.*[bw*Q*Sw;cbar*Q*Sw;bw*Q*Sw];

%-------------------------8. Engine Force & Moment-------------------------

%Assuming engine thrust is aligned with body frame
% F1=u4*m*g;
% F2=u5*m*g;
% FE1_b = [F1;
%         0;
%         0];
% FE2_b = [F2;
%         0;
%         0];

% FE_b = FE1_b + FE2_b;

%propeller thrust - Fitzpatrick model for ACP 25x12.5E
km = 37.42;% motor constant [m/s] torque/sqrt(power) - ne znam
Cp = 0.57;% efficientcy factor
D = 0.635; % diameter [m]
Sp = pi*(D^2/4); % disc area
d_t1 = u4/Tp; % range [0,1]
Vd1 = Va + d_t1*(km-Va);
Tp1 = 0.5*rho*Sp*Cp*Vd1*(Vd1-Va);
d_t2 = u5/Tp; % range [0,1]
Vd2 = Va + d_t2*(km-Va);
Tp2 = 0.5*rho*Sp*Cp*Vd2*(Vd2-Va);
FE1_b = [Tp1*cos(-5*pi/180);
        0;
        sin(5*pi/180)];
FE2_b = [Tp2*cos(-5*pi/180);
        0;
        sin(5*pi/180)];
FE_b = FE1_b + FE2_b;
%Now engine moment due to offset of thrust from CG
mew1 = [Xcg-Xapt1;
        Yapt1-Ycg;
        Zcg-Zapt1];

mew2 = [Xcg-Xapt2;
        Yapt2-Ycg;
        Zcg-Zapt2];

MEcg1_b = cross(mew1,FE1_b);
MEcg2_b = cross(mew2,FE2_b);

MEcg_b = MEcg1_b + MEcg2_b;
%----------------------------9. Gravity Effect-----------------------------
%Calculate gravitational forces in body frame. This causes no moment
%about--cg

g_b = [-g*sin(x8);g*cos(x8)*sin(x7);g*cos(x8)*cos(x7)];

Fg_b = m*g_b;

%----------------------------10. State Derivatives------------------------
%Inertia Matrix
Ib = [39.71 0 8.97;
      0 85.51 0;
      8.97 0 114.39];
invIb = inv(Ib);
           
F_b = Fg_b + FE_b + FA_b;
x1to3dot = (1/m)*F_b - cross(wbe_b,V_b);

Mcg_b = MAcg_b + MEcg_b;
x4to6dot = invIb*(Mcg_b - cross(wbe_b,Ib*wbe_b));

% H_phi = [1 sin(x7)*tan(x8) cos(x7)*tan(x8);
%          0 cos(x7) -sin(x7);
%          0 sin(x7)/cos(x8) cos(x7)/cos(x8)];
%MFS
H_phi = [1 0 0;
         0 cos(alpha) -sin(beta);
         0 sin(alpha) cos(alpha)];
x7to9dot = H_phi*wbe_b;
x10to12dot = rot_body_to_ned(X(1:9));
XDOT = [x1to3dot;
    x4to6dot;
    x7to9dot;
    x10to12dot;
    LD_ratio;
    F_b;
    Mcg_b;
%     CL;
    CQ;
%     CD;
    Cl;
    Cm;
    Cn;
    alpha;
    beta;
    CL_w_OGE;
    CL_h_OGE;
    CL_w_IGE;
    CL_h_IGE;
    CD_iw_IGE;
    CD_ih_IGE;
    Fg_b;
    FE_b; 
    FA_b];


