function double_pendulum_minimal
% ------------------------------------------------------------
%   4-state  ODE
% ------------------------------------------------------------

%%  -------------------------------------------------------
p.m    = 1;                  % kg
p.L    = 1;                  % m
p.J    = p.m * p.L^2/12;     % kg·m²
p.g    = 9.81;               % m/s²   
p.b    = 10;                 % N·m·s
p.uamp = 5;                  % N·m    
p.ufrq = 5;                  % rad/s

%% initial x = [θ₁; θ₂(abs); ω₁; ω₂(abs)] -----------------------
x0 = [0; 0; 0; 0];           

%% time ---------------------------------------------------
tEnd = 10;
opt  = odeset('RelTol',1e-8,'AbsTol',1e-11);
[t,X] = ode15s(@(t,x) rhs_minimal(t,x,p), [0 tEnd], x0, opt);

%% ——  —— -------------------------
th1 = X(:,1);                     % ab θ₁
th2 = X(:,2) - X(:,1);            % real θ₂   
w1  = X(:,3);                     
w2  = X(:,4) - X(:,3);            % real ω₂  

figure('Color','w');
plot(t, th1,'b',  t, th2,'c--', ...
     t, w1 ,'r',  t, w2 ,'m--','LineWidth',1.3);
grid on
xlabel('time  [s]')
ylabel('\theta  [rad]  /  \omega  [rad·s^{-1}]')
legend('\theta_1','\theta_2  (rel)', ...
       '\omega_1','\omega_2 (rel)','Location','best')
title('Double pendulum – minimal coords  (matches maximal / Simscape)')

assignin('base','t_min',t);
assignin('base','X_min',X);
end
% ========================================================================
function dx = rhs_minimal(t,x,p)
% 4-state RHS  ——  x = [θ₁; θ₂(abs); ω₁; ω₂(abs)]
th1 = x(1);  th2 = x(2);
w1  = x(3);  w2 = x(4);

m = p.m;  L = p.L;  J = p.J;  g = p.g;  b = p.b;
u = p.uamp * sin(p.ufrq*t);          %input

c12 = cos(th1 - th2);

% --------  M(θ) --------------------------------------
M11 = J + 5/4 * m*L^2;
M12 = 1/2 * m*L^2 * c12;
M22 = J + 1/4 * m*L^2;
detM = M11*M22 - M12^2;             

% --------G --------------------------------
h1 = 3/2 * m*g*L * cos(th1);
h2 = 1/2 * m*g*L * cos(th2);

% -------- damping ---------------------------------------
tau1 = u + b * (w2 - 2*w1);         % 1
tau2 =     b * (w1 -   w2);         % 2

rhs1 = tau1 - h1;
rhs2 = tau2 - h2;

% --------  α = M⁻¹ (τ − h) ----------------------------
a1 = ( M22*rhs1 - M12*rhs2) / detM;
a2 = (-M12*rhs1 + M11*rhs2) / detM;

dx = [ w1 ;
       w2 ;
       a1 ;
       a2 ];
end
