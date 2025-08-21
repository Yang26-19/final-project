function [tlog, th1_log, rel12_log, rel23_log] = simulate_method_full_3link(alpha_speed, T)
% simulate_method_full_3link
% 3-link serial chain: world ↔ link1 ↔ link2 ↔ link3
% Unknowns X = [x1dd y1dd th1dd x2dd y2dd th2dd x3dd y3dd th3dd  λ1 λ2 λ3 λ4 λ5 λ6]^T
% λ1,λ2: base (world↔1);  λ3,λ4: joint 1 (1↔2);  λ5,λ6: joint 2 (2↔3)
% Outputs: th1 (absolute), (th2 - th1), (th3 - th2)
%
% Note on constraints:
% - Base constraints (world↔link1) are absolute: the world point has zero acceleration,
%   so the residual equals the endpoint acceleration of link1 (plus Adot*dq).
% - Joint constraints are relative: residual is "left endpoint minus right endpoint"
%   (plus Adot*dq). This is why base residuals look sign-different from joint ones.

if nargin<1, alpha_speed = 3000; end
if nargin<2, T = 4; end

% ---- Parameters ----
g0 = 9.81;   % gravity
b  = 0.5;    % joint viscous damping (relative)
L  = 1.0;    % link length
m  = 1.0;    % link mass
J  = L^2/12; % link moment of inertia about COM (slender rod)
Lh = 0.5*L;  % half-length to endpoints

% ---- Time ----
dt    = 1e-4;
steps = ceil(T/dt);
tlog  = (0:steps-1)*dt;
Twait = 0.02;  % wait before integrating outer angles (lets inner loop settle)

% ---- Outer (angle) states ----
th1=0; th2=0; th3=0;     % angles
w1=0;  w2=0;  w3=0;      % angular rates

% ---- Inner unknowns (accelerations and lambdas) ----
X = zeros(15,1);

% ---- Logs ----
th1_log   = zeros(1,steps);
rel12_log = zeros(1,steps);
rel23_log = zeros(1,steps);

for k = 1:steps
    % Trig
    s1=sin(th1); c1=cos(th1);
    s2=sin(th2); c2=cos(th2);
    s3=sin(th3); c3=cos(th3);

    % Unpack X
    x1dd=X(1);  y1dd=X(2);  th1dd=X(3);
    x2dd=X(4);  y2dd=X(5);  th2dd=X(6);
    x3dd=X(7);  y3dd=X(8);  th3dd=X(9);
    lam1=X(10); lam2=X(11); lam3=X(12); lam4=X(13); lam5=X(14); lam6=X(15);

    % ---- Relative viscous torques τi at joints (w0=0, w4=0) ----
    % Each tau is the net damping torque on the link due to relative rates.
    tau1 = -b*w1            + b*(w2 - w1);
    tau2 = -b*(w2 - w1)     + b*(w3 - w2);
    tau3 = -b*(w3 - w2);

    % ---- Adot*dq terms (Coriolis/centripetal-like) ----
    % Base constraint is absolute (to world, taken as inertial), joints are relative.
    ad1 = +Lh*c1*w1^2;                        % base x
    ad2 = +Lh*s1*w1^2;                        % base y
    ad3 = -Lh*(c1*w1^2 + c2*w2^2);            % joint1 x (left minus right)
    ad4 = -Lh*(s1*w1^2 + s2*w2^2);            % joint1 y
    ad5 = -Lh*(c2*w2^2 + c3*w3^2);            % joint2 x
    ad6 = -Lh*(s2*w2^2 + s3*w3^2);            % joint2 y

    % ================= "Top" equations: e1..e9 (link dynamics + gravity + joint/base forces) =================
    % link1
    e1 = m*x1dd - (lam1 + lam3);
    e2 = m*y1dd - (lam2 + lam4) + m*g0;
    e3 = J*th1dd - ( Lh*s1*lam1 - Lh*c1*lam2 - Lh*s1*lam3 + Lh*c1*lam4 ) - tau1;

    % link2 (left +lam3,lam4; right -lam5,lam6)
    e4 = m*x2dd + lam3 - lam5;
    e5 = m*y2dd + lam4 - lam6 + m*g0;
    e6 = J*th2dd - ( -Lh*s2*lam3 + Lh*c2*lam4 ) ...
                  - ( -Lh*s2*lam5 + Lh*c2*lam6 ) - tau2;

    % link3 (only the left joint2 acts)
    e7 = m*x3dd + lam5;
    e8 = m*y3dd + lam6 + m*g0;
    e9 = J*th3dd - ( -Lh*s3*lam5 + Lh*c3*lam6 ) - tau3;

    % ================= "Bottom" constraint residuals: c1r..c6r =================
    % base (world↔1): absolute constraint (world acceleration = 0)
    c1r = x1dd + Lh*s1*th1dd + ad1;
    c2r = y1dd - Lh*c1*th1dd + ad2;

    % joint1 (1↔2): relative (left minus right)
    c3r = (x1dd - Lh*s1*th1dd) - (x2dd + Lh*s2*th2dd) + ad3;
    c4r = (y1dd + Lh*c1*th1dd) - (y2dd - Lh*c2*th2dd) + ad4;

    % joint2 (2↔3): relative (left minus right)
    c5r = (x2dd - Lh*s2*th2dd) - (x3dd + Lh*s3*th3dd) + ad5;
    c6r = (y2dd + Lh*c2*th2dd) - (y3dd - Lh*c3*th3dd) + ad6;

    % ================= Gradient g = K^T r : g1..g15 =================
    % ddq columns
    g1 =  m*e1 + c1r + c3r;                                           % x1dd
    g2 =  m*e2 + c2r + c4r;                                           % y1dd
    g3 =  J*e3 + ( +Lh*s1)*c1r + ( -Lh*c1)*c2r ...
              + ( -Lh*s1)*c3r + ( +Lh*c1)*c4r;                        % th1dd

    g4 =  m*e4 - c3r + c5r;                                           % x2dd
    g5 =  m*e5 - c4r + c6r;                                           % y2dd
    g6 =  J*e6 + ( -Lh*s2)*(c3r + c5r) + ( +Lh*c2)*(c4r + c6r);       % th2dd

    g7 =  m*e7 - c5r;                                                 % x3dd
    g8 =  m*e8 - c6r;                                                 % y3dd
    g9 =  J*e9 + ( -Lh*s3)*c5r + ( +Lh*c3)*c6r;                       % th3dd

    % lambda columns (using only the Top rows)
    g10 = -( e1 + Lh*s1*e3 );                                         % λ1 (base x)
    g11 = -( e2 - Lh*c1*e3 );                                         % λ2 (base y)
    g12 = -( e1 - Lh*s1*e3 - e4 - Lh*s2*e6 );                         % λ3 (joint1 x)
    g13 = -( e2 + Lh*c1*e3 - e5 + Lh*c2*e6 );                         % λ4 (joint1 y)
    g14 = -( e4 - Lh*s2*e6 - e7 - Lh*s3*e9 );                         % λ5 (joint2 x)
    g15 = -( e5 + Lh*c2*e6 - e8 + Lh*c3*e9 );                         % λ6 (joint2 y)

    gvec = [g1;g2;g3; g4;g5;g6; g7;g8;g9; g10;g11;g12;g13;g14;g15];

    % ---- Inner gradient descent step ----
    X = X - dt*alpha_speed*gvec;

    % ---- Outer double integration (after Twait) ----
    if tlog(k) >= Twait
        th1dd = X(3);  th2dd = X(6);  th3dd = X(9);
        w1 = w1 + dt*th1dd;  w2 = w2 + dt*th2dd;  w3 = w3 + dt*th3dd;
        th1 = th1 + dt*w1;   th2 = th2 + dt*w2;   th3 = th3 + dt*w3;
    end

    % ---- Logs ----
    th1_log(k)   = th1;
    rel12_log(k) = th2 - th1;
    rel23_log(k) = th3 - th2;
end

% ---- Plot: θ1 and relative angles ----
figure('Color','w','Name','3-link: th1 & relatives');
plot(tlog, th1_log,'LineWidth',1.35); hold on;
plot(tlog, rel12_log,'--','LineWidth',1.35);
plot(tlog, rel23_log,':','LineWidth',1.35);
grid on;
ylabel('\theta (rad)');
xlabel('t (s)');
title(sprintf('\\alpha=%.0f, T=%.1fs, T_{wait}=%.3g s',alpha_speed,T,Twait));
legend('\theta_1','\theta_2-\theta_1','\theta_3-\theta_2','Location','best');

% ---- Absolute angles (recover th2, th3) ----
th2_log = th1_log + rel12_log;
th3_log = th2_log + rel23_log;

figure('Color','w','Name','3-link: absolute thetas');
plot(tlog, th1_log,'LineWidth',1.35); hold on;
plot(tlog, th2_log,'--','LineWidth',1.35);
plot(tlog, th3_log,':','LineWidth',1.35);
grid on;
ylabel('\theta (rad)');
xlabel('t (s)');
legend('\theta_1','\theta_2','\theta_3','Location','best');
title('Absolute angles');
end
