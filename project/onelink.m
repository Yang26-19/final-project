function simulate_gradient_with_outer(alpha_speed, Twait)

    % ---- Physical Parameters ----
    g0 = 9.81; 
    b  = 0.5; 
    I  = 1/12;

    % ---- Initial State ----
    X = zeros(5,1);        % [xdd; ydd; thdd; lam1; lam2]
    theta = 0;         % 15 degrees initial angle
    dtheta = 0;            % angular velocity

    % ---- Time Settings ----
    dt = 1e-4;
    steps = ceil(Twait/dt);
    tlog = (0:steps-1)*dt;

    % ---- Logging ----
    X_log = zeros(5, steps);
    res_log = zeros(1, steps);
    theta_log = zeros(1, steps);
    dtheta_log = zeros(1, steps);

    % ---- Gradient Descent + Outer Dynamics ----
    for k = 1:steps
        t = tlog(k);
        s = sin(theta); 
        c = cos(theta);

        % External forces
        f1 = 0;
        f2 = -g0;
        f3 = -b*dtheta;
        f4 = -0.5*dtheta^2 * c;
        f5 = -0.5*dtheta^2 * s;

        % Current values
        xdd  = X(1); 
        ydd  = X(2); 
        thdd = X(3); 
        lam1 = X(4); 
        lam2 = X(5);

        % Error vector: e = A*X - b
        e1 =  xdd            + lam1                          - f1;
        e2 =  ydd                        + lam2              - f2;
        e3 =          I*thdd + 0.5*s*lam1 - 0.5*c*lam2       - f3;
        e4 =  xdd + 0.5*s*thdd                               - f4;
        e5 =  ydd - 0.5*c*thdd                               - f5;

        % Gradient: g = A^T * e
        g1 = e1 + e4;
        g2 = e2 + e5;
        g3 = I*e3 + 0.5*s*e4 - 0.5*c*e5;
        g4 = e1 + 0.5*s*e3;
        g5 = e2 - 0.5*c*e3;
        gvec = [g1; g2; g3; g4; g5];

        % Gradient descent update (inner loop)
        X = X - dt * alpha_speed * gvec;

        % Outer dynamics update
        if t >= 0.02
    thdd  = X(3);
    dtheta = dtheta + dt * thdd;
    theta  = theta  + dt * dtheta;
        end


        % Logging
        X_log(:, k) = X;
        res_log(k) = norm([e1; e2; e3; e4; e5]);
        theta_log(k) = theta;
        dtheta_log(k) = dtheta;
    end

    % ---- Plotting Results ----
    figure;

    subplot(4,1,1);
    plot(tlog, theta_log, 'LineWidth', 1.5); grid on;
    ylabel('\theta (rad)');
    title(sprintf('Full Simulation with Outer Dynamics (\\alpha_{speed}=%.0f)', alpha_speed));

    subplot(4,1,2);
    plot(tlog, dtheta_log, 'LineWidth', 1.5); grid on;
    ylabel('d\\theta/dt (rad/s)');

    subplot(4,1,3);
    plot(tlog, X_log(1,:), 'r', tlog, X_log(2,:), 'g', tlog, X_log(3,:), 'b', 'LineWidth', 1.5);
    legend('xdd', 'ydd', 'thdd'); grid on;
    ylabel('ddq');

    subplot(4,1,4);
    semilogy(tlog, res_log, 'k', 'LineWidth', 1.2); grid on;
    xlabel('Time (s)');
    ylabel('||A X - b||');
end

% Run the simulation
simulate_gradient_with_outer(3000, 10);
