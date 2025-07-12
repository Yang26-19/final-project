function link1
    
    m = 1.0;
    L = 1.0;
    g = 10.0;
    b = 0.1;
    I_COM = (1/3) * m * L^2;

    % initial ：theta = 0
    theta0 = 0;
    omega0 = 0;
    x0 = (L/2) * cos(theta0);
    y0 = (L/2) * sin(theta0);
    dx0 = - (L/2) * sin(theta0) * omega0;
    dy0 = (L/2) * cos(theta0) * omega0;

    z0 = [x0; y0; theta0; dx0; dy0; omega0];

    tspan = [0 10];

    [t, Z] = ode45(@(t, z) odeFunc(t, z, m, L, g, b, I_COM), tspan, z0);

    x = Z(:,1); y = Z(:,2); theta = Z(:,3);
    dx = Z(:,4); dy = Z(:,5); omega = Z(:,6);

    % Plot
    figure;
    subplot(3,1,1);
    plot(t, theta, 'LineWidth', 1.5);
    ylabel('\theta (rad)');
    title('\theta vs Time'); grid on;

    subplot(3,1,2);
    plot(t, omega, 'LineWidth', 1.5);
    ylabel('\omega (rad/s)');
    title('Angular Velocity'); grid on;

    subplot(3,1,3);
    plot(t, x, 'b', t, y, 'r', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Position (m)');
    legend('x', 'y'); title('Center of Mass Position'); grid on;
end

function dz = odeFunc(~, z, m, L, g, b, I_COM)
    x = z(1); y = z(2); theta = z(3);
    dx = z(4); dy = z(5); omega = z(6);

    % Step 4: Angular acceleration from dynamics
   ddtheta = -(m * g * (L/2) * cos(theta) + b * omega) / I_COM;


    % Step 5: New geometry-based expressions for COM acceleration
    ddx = - (L/2) * (cos(theta) * omega^2 + sin(theta) * ddtheta);
    ddy = (L/2) * (cos(theta) * ddtheta - sin(theta) * omega^2);

    dz = zeros(6,1);
    dz(1) = dx;
    dz(2) = dy;
    dz(3) = omega;
    dz(4) = ddx;
    dz(5) = ddy;
    dz(6) = ddtheta;
end
