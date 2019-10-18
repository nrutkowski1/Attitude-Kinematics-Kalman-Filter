function kalman_filter_attitude_kinematics()

    load sensorlog_20190923_151441.mat...
        Orientation AngularVelocity 
    
%     load attitude_kinematics_synthetic_data1.mat...
%         euler_angles_0 euler_angles_noisy euler_angles_true...
%         omega_iee_noisy omega_iee_true...
%         std_dev_gyro std_dev_mag time_pts
    TT = synchronize(Orientation, AngularVelocity, 'regular', 'spline','TimeStep', seconds(0.01));
    
    orientation = [TT.Z_Orientation*pi/180 TT.Y_Orientation*pi/180 TT.X_Orientation*pi/180]';
    rate_gyros = [TT.Z_AngularVelocity TT.Y_AngularVelocity TT.X_AngularVelocity]';
    
    time_pts = TT.Timestamp;
    
%     dt = time_pts(2) - time_pts(1);
    
%     orientation = euler_angles_noisy;
%     rate_gyros = omega_iee_noisy; 
    var_rgy = 0.2932;
    var_mag = 2.385*pi/180;
%     var_rgy = std_dev_gyro^2;
%     var_mag = std_dev_mag^2;
    % changes based on problem 
    n_states = 3;
    C = eye(n_states);
    % only here if constant
    R = var_mag*eye(3);
    G = var_rgy*eye(3);
    
        
    % initialize 
    psi0 = orientation(1, 1);
    theta0 = orientation(2, 1);
    phi0 = orientation(3, 1);
    
    omega_aee_t = rate_gyros(:, 1);
    
    x_hat = [psi0; theta0; phi0];
%     x_hat = euler_angles_0;
    P = R; 
    dt = 0.01;
    % plotting
    x_hat_rec = zeros(3, numel(time_pts));
    x_hat_rec(:, 1) = x_hat; 
    
    for n = 2:numel(time_pts)

        % A and B Matrics
        psi = orientation(1, n);
        theta = orientation(2, n);
        phi = orientation(3, n);

        A = eye(n_states);
        H321 = [0 sin(phi) cos(phi); 
                0 cos(phi)*cos(theta) -sin(phi)*cos(theta);
                cos(theta) sin(phi)*sin(theta) cos(phi)*sin(theta)];
            
        B = (dt/cos(theta))*H321;
            
        Q = B*G*B';
        
        % predictive update
        omega_aee_t = rate_gyros(:, n);
        
        x_minus = rk4_step(x_hat);
        
        P_minus = A*P*A' + Q;
        
        % Kalman gain
        K = P_minus*C' / (C*P_minus*C' + R);
        
        % measurement update
        z = [psi; theta; phi];
        x_hat = x_minus + K*(z - C*x_minus);
        P = (eye(n_states) - K*C)*P_minus; 
        
        % sum of diagonal elements 
        trace(P);
        
        x_hat_rec(:, n) = x_hat;
    end
    
    figure(1);
    subplot(311);
    plot(time_pts, x_hat_rec(1, :), 'b-', time_pts, orientation(1, :), 'r-');
    ylabel('\psi (rad)');
    legend('calculated', 'measured')
    subplot(312);
    plot(time_pts, x_hat_rec(2, :), 'b-', time_pts, orientation(2, :), 'r-');
    ylabel('\theta (rad)')
    subplot(313);
    plot(time_pts, x_hat_rec(3, :), 'b-', time_pts, orientation(3, :), 'r-');
    ylabel('\phi (rad)')
    xlabel('Time (s)')   
    
    function x_dot = f_x(x_t)
        
        e_thta = x_t(2);
        e_phi = x_t(3);
       
        x_dot = [-sin(e_thta) 0 1; 
                sin(e_phi)*cos(e_thta) cos(e_phi) 0;
                cos(e_phi)*cos(e_thta) -sin(e_phi) 0] \ omega_aee_t;
        
    end

    function x_t_plus_dt = rk4_step(x_t)
       
        k1 = dt * f_x(x_t);
        k2 = dt * f_x(x_t + 0.5*k1);
        k3 = dt * f_x(x_t + 0.5*k2);
        k4 = dt * f_x(x_t + k3);
        
        x_t_plus_dt = x_t + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4;       
    end

end