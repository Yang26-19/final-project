function double_pendulum_maximal
% ------------------------------------------------------------
%  Planar 2-link pendulum (maximal coords, matches Simscape output)
% ------------------------------------------------------------

%% -------- parameters I use --------------------------
par.L = 1;           par.m = 1;          par.J = par.L^2/12;
par.g = 9.81;        par.b = 10;         % viscous coef (both joints)
par.uamp = 5;        par.ufrq = 5;       % drive on joint-1
tEnd = 10;           wantDbg = false;

%% -------- initial pose & velocities -------------------------
th1_0 = 0;                        % both links horizontal
th2rel_0 = 0;
th2abs_0 = th1_0 + th2rel_0;

L = par.L;
q0 = [ L/2*cos(th1_0);            % x1
       L/2*sin(th1_0);            % y1
       th1_0;                     % θ1 abs
       L*cos(th1_0)+L/2*cos(th2abs_0); % x2
       L*sin(th1_0)+L/2*sin(th2abs_0); % y2
       th2abs_0 ];                % θ2 abs
dq0 = zeros(6,1);
y0  = [q0; dq0];

opts = odeset('RelTol',1e-8,'AbsTol',1e-11);
[t,Y] = ode15s(@(t,y) rhs(t,y,par), [0 tEnd], y0, opts);

%% -------- get joint variables --------------------------------
th1 = Y(:,3);                    % θ1 abs
th2abs = Y(:,6);                 % θ2 abs
th2 = th2abs - th1;              % θ2 relative

w1 = Y(:,9);                     % ω1 abs
w2abs = Y(:,12);                 % ω2 abs
w2 = w2abs - w1;                 % ω2 relative

%% -------- plot -----------------------------------------------
figure('Color','w');
plot(t, th1,'b', t, th2,'c--', t, w1,'r', t, w2,'m--','LineWidth',1.2);
legend('\theta_1 (rad)', '\theta_2 (rad)', ...
       '\omega_1 (rad/s)', '\omega_2 (rad/s)', 'Location','best');
xlabel('time  [s]');
ylabel('value  [rad / rad·s^{-1}]');
title('Double pendulum – maximal coords (matches Simscape)');
grid on;

%% optional ----------------------------------------------------
if wantDbg
    [h,lambda] = constraint_and_lambda(Y,par);
    fprintf('max |h_i| = %.3e\n', max(abs(h),[],'all'));
end
end
% =================================================================
function dy = rhs(t,y,p)
q  = y(1:6);  dq = y(7:12);

% shorthand ---------------------------------------------------
th1 = q(3);  th2 = q(6);
s1 = sin(th1); c1 = cos(th1);
s2 = sin(th2); c2 = cos(th2);
L  = p.L;

% Jacobian A --------------------------------------------------
A = [ 1 0  0.5*L*s1   0 0 0;
      0 1 -0.5*L*c1   0 0 0;
      1 0 -0.5*L*s1  -1 0 -0.5*L*s2;
      0 1  0.5*L*c1   0 -1  0.5*L*c2];

% mass, gravity ----------------------------------------------
M = diag([p.m p.m p.J p.m p.m p.J]);
G = [0; p.m*p.g; 0;   0; p.m*p.g; 0];   % Y-down positive

% --- joint torques ------------------------------------------
w1      = dq(3);
w2abs   = dq(6);
w2rel   = w2abs - w1;

tau1_damp = -p.b*w1 + p.b*w2rel;   % from both joints
tau2_damp = -p.b*w2rel;

u = p.uamp * sin(p.ufrq*t);        % drive on joint-1

Q = [0;0;  u + tau1_damp;  0;0;  tau2_damp];

% --- mixed linear system ------------------------------------
K   = [M  -A';
       A   zeros(4)];
rhs = [Q - G;
       -A*dq];                     % first-order Baumgarte γ=−A·dq

z   = K\rhs;
ddq = z(1:6);

dy = [dq ; ddq];
end
% =================================================================
function [h,lambda] = constraint_and_lambda(Y,p)
q   = Y(:,1:6);
th1 = q(:,3);   th2 = q(:,6);
s1=sin(th1); c1=cos(th1);
s2=sin(th2); c2=cos(th2);
L = p.L;

h = [ q(:,1) - 0.5*L.*c1, ...
      q(:,2) - 0.5*L.*s1, ...
      q(:,1)+0.5*L.*c1 - q(:,4)+0.5*L.*c2, ...
      q(:,2)+0.5*L.*s1 - q(:,5)+0.5*L.*s2 ];

lambda = [];         % 只看误差时可以不用求 λ
end
