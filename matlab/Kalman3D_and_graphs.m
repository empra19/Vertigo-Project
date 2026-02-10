
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gpsdata = evalin('base', 'gpsdata');
imudata = evalin('base', 'imudata');
accel_ned = evalin('base', 'accel_ned');

% Convert GPS data into metres from starting position
gpstime = gpsdata(:,1);
[east, north] = ll2utm(gpsdata(:,4), gpsdata(:,3));
dalt = -(gpsdata(:,5) - gpsdata(1,5)); % Down is positive
deast = east - east(1);
dnorth = north - north(1);

% Interpolate GPS data onto the IMU data
time = imudata(:,1);
deast = interp1(gpstime, deast, time, 'pchip', 'extrap');
dnorth = interp1(gpstime, dnorth, time, 'pchip', 'extrap');
dalt = interp1(gpstime, dalt, time, 'pchip', 'extrap');

% Run Kalman filter
Nstep = length(time);

% State vector: [s_n, s_e, s_d, v_n, v_e, v_d, a_n, a_e, a_d]
x = zeros(9, 1);
P = 1e-3 * eye(9);

dt = time(2) - time(1); % This is still a really bad plan

% Process noise
Q = diag([1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-3 1e-3 1e-3]);
% Measurement noise
gps_var = 1e-6;
imu_var = 1e-3;
R = diag([gps_var gps_var gps_var*10 imu_var imu_var imu_var]);

% Discrete time model
% s_n(k+1) = s_n(k) + v_n(k) * dt
% v_n(k+1) = v_n(k) + a_n(k) * dt
% a_n(k+1) = a_n(k)
F = [1 0 0 dt 0  0  0  0  0; ...
     0 1 0 0  dt 0  0  0  0; ...
     0 0 1 0  0  dt 0  0  0; ...
     0 0 0 1  0  0  dt 0  0; ...
     0 0 0 0  1  0  0  dt 0; ...
     0 0 0 0  0  1  0  0  dt; ...
     0 0 0 0  0  0  1  0  0; ...
     0 0 0 0  0  0  0  1  0; ...
     0 0 0 0  0  0  0  0  1];
 
H = [1 0 0 0 0 0 0 0 0; ...
     0 1 0 0 0 0 0 0 0; ...
     0 0 1 0 0 0 0 0 0; ...
     0 0 0 0 0 0 1 0 0; ...
     0 0 0 0 0 0 0 1 0; ...
     0 0 0 0 0 0 0 0 1];

%  Now in vtg_load_data_and_transform
% % Remove gravity from NED accel, leaving linear accels in NED frame
% accel_ned(:,3) = accel_ned(:,3) - 1;
% accel_ned = accel_ned .* 9.81;

% Sensors
sensors = @(i) [dnorth(i);
                deast(i);
                dalt(i);
                accel_ned(i, 1:3)'];
                
kal_x_stor = zeros(9, Nstep);

h = waitbar(0, 'Running Kalman filter...');
for i = 1:Nstep
    % Predict
    xp = F * x; % no inputs
    Pp = F * P * F' + Q;
    
    % Update
    y = sensors(i) - H * x;
    S = H * P * H' + R;
    K = Pp * H' * inv(S);
    x = xp + K * y;
    P = (eye(9) - K * H) * Pp;
    
    % Store
    kal_x_stor(:, i) = x;
    
    % Waitbar
    if mod(i/Nstep, 0.1) == 0
        waitbar(i/Nstep, h);
    end
end
close(h);

% Plot results
clf;
% subplot(2,1,1);
%plot3(kal_x_stor(2, :), kal_x_stor(1, :), -kal_x_stor(3, :));
%hold on
plot3(deast, dnorth, -dalt);
xlabel('North (m)');
ylabel('East (m)');
%legend ('Fusion Position', 'GPS')
grid on
zlabel ('Altitude (m)');
%view(90,90);

figure
% subplot(2,1,2);
plot(time, -kal_x_stor(3, :));
hold on
plot(time, -dalt);
xlabel('Time (s)');
ylabel('Altitude (m)')





%end % function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Yaw angle
Phi_radians = atan2(2.*((quatdata(:,3).*quatdata(:,6))+(quatdata(:,4).*quatdata(:,5))),1 - 2.*(((quatdata(:,3).^2)+(quatdata(:,4).^2))));
Phi = (Phi_radians)*180./pi;
Yaw = Phi;
%roll angle
Theta_radians = asin (2.*(quatdata(:,4).*quatdata(:,6) -   quatdata(:,3).*quatdata(:,5)));
Theta = (Theta_radians).*180./pi;
Roll = Theta;
% Pitch angle
Psi_radians = atan2(2.*((quatdata(:,5).*quatdata(:,6))+(quatdata(:,3).*quatdata(:,4))),1 - 2.*(((quatdata(:,4).^2)+(quatdata(:,5).^2))));
Psi = (Psi_radians).*180./pi;
Pitch = Psi;
quiver_length = 3;
[rollx,rolly] = pol2cart(((Theta_radians)),quiver_length);
[pitchx,pitchy] = pol2cart(((Psi_radians+pi)),quiver_length);
[yawx,yawy] = pol2cart((Phi_radians),quiver_length);

%Cartesian vectors for roll, pitch and yaw with time
roll_pitch_yaw_t = [rollx rolly pitchx pitchy yawx yawy quatdata(:,1)];
Angles = [Roll Pitch Yaw quatdata(:,1)];

prompt = 'What time do you wish to start analysis from? ';
window_start = input(prompt);
prompt = 'What time do you wish to end the analysis? ';
window_end = input(prompt);

%Total time samples
total_time = [imudata(:,1) ; gpstime];
%sort in time chronological order
stotal_time = sort(total_time);






% Extract the bit of data we want to look at
tstartidx = find(imudata(:,1) > window_start, 1);
tendidx = find(imudata(:,1) > window_end, 1);
Time_window = time(tstartidx:tendidx,1);

tstartquat = find(quatdata(:,1) > window_start, 1);
tendquat = find(quatdata(:,1) > window_end, 1);



kal_x_storT =kal_x_stor';
smooth_north_position = smooth (kal_x_storT(tstartidx:tendidx,1));
smooth_east_position = smooth(kal_x_storT(tstartidx:tendidx,2));
smooth_down_position = smooth(kal_x_storT(tstartidx:tendidx,3));
north_vel_elements = smooth (kal_x_storT(tstartidx:tendidx,4));
east_vel_elements = smooth (kal_x_storT(tstartidx:tendidx,5));
down_vel_elements = smooth (kal_x_storT(tstartidx:tendidx,6));
north_accel_elements = smooth(kal_x_storT(tstartidx:tendidx,7));
east_accel_elements = smooth(kal_x_storT(tstartidx:tendidx,8));
down_accel_elements = smooth (kal_x_storT(tstartidx:tendidx,9));
roll_pitch_yaw_twindow = roll_pitch_yaw_t (tstartquat:tendquat,:);
Angles_window = Angles(tstartquat:tendquat,:);
% quat_data = quatdata(tstartidx:tendidx, :);
% euldata_window = euldata (tstartidx:tendidx, :);

NEDT = [smooth_north_position smooth_east_position smooth_down_position Time_window];

% Find the spot velocities and accelerations as taken from the kalman position
VelN =  smooth((movmean (diff(NEDT(:,1))./0.01,100)));
VelE = smooth( (movmean (diff(NEDT(:,2))./0.01,100)));
VelD = smooth((movmean (diff(NEDT(:,3))./0.01,100)));
accelN = smooth(movmean (diff(VelN)./0.01,100));
accelE =smooth( movmean (diff(VelE)./0.01,100));
accelD = smooth(movmean (diff(VelD)./0.01,100));

NEDT_vels = [VelN VelE VelD Time_window(1:end-1,:)];

NEDT_accels = [accelN accelE accelD Time_window(1:end-2,:)];

%calculate the magnitude of the velocity and acceleration
Abs_accel = sqrt((accelN.^2) + (accelE.^2) + (accelD.^2));
Abs_vel = sqrt ((VelN.^2) + (VelE.^2) + (VelD.^2));
Abs_pos = sqrt ((smooth_north_position.^2) + (smooth_east_position.^2) + (smooth_north_position.^2));

z = 1;
i = 1;
% Merging angles with time 
%Angles_time has N E D T R P Y

Angles_time = [];
    
       Angles_time = horzcat(NEDT(:,1:3), Angles_window);

%joining all data position and roll pitch yaw and equalising matrix sizes
all_data = [];
    for a = 1 :length(NEDT)
        all_data (a,:) = horzcat(NEDT (a,1:3), roll_pitch_yaw_twindow(a,:));
    end
    

z = 1;
%all_data holds [N E D Rx Ry Px Py Yx Yy T]

universal_time = all_data(:,10);
UT = universal_time;
z = 1;
%join all data velocity and roll pitch yaw and equalise matrix size
for z = 1 :length(NEDT_vels)    
    for  i = 1: length(roll_pitch_yaw_t(:,1))

    
    if roll_pitch_yaw_t(i,7) == NEDT_vels(z,4);
        all_data_vel (i,:) = horzcat(NEDT_vels (z,1:3), roll_pitch_yaw_t(i,:));
    end
end
end


% joining all position data accelerations data
for z = 1 :length(NEDT(1:end-2,:))    
    for  i = 1: length(UT(1:end-2,:))

    
    if NEDT_accels (i,4) == NEDT_vels(z,4);
        all_data_and_accels (i,:) = horzcat(NEDT_vels (z,1:3), NEDT_accels(i,:));
    end
end
end


% GRAPHS GRAPHS GRAPHS
% GRAPHS GRAPHS GRAPHS
% GRAPHS GRAPHS GRAPHS
% GRAPHS GRAPHS GRAPHS
% GRAPHS GRAPHS GRAPHS
z = 1;

%decimate_rate  = dr  
% You can reduce the number of data points by increasing this value.
dr = 10;

figure;
plot (Angles_time(:,7), Angles_time(:,6), Angles_time(:,7), Angles_time(:,4),Angles_time(:,7), Angles_time(:,5));
xlabel('Time (s)');
ylabel('Angle/degrees');
legend('Yaw', 'Roll', 'Pitch');

figure;
subplot (3,1,1)

plot(Angles_time(:,7), Angles_time(:,4));
xlabel('Time (s)');
ylabel('Angle/ deg');
legend('Roll NED');

subplot (3,1,2)

plot(Angles_time(:,7), Angles_time(:,5));
xlabel('Time (s)');
ylabel('Angle/ deg');
legend('Pitch NED');

subplot (3,1,3)

plot(Angles_time(:,7), Angles_time(:,6));
xlabel('Time (s)');
ylabel('Angle/ deg');
legend('Yaw NED');


%quiver plot yaw

% 
% %quiver plot roll

figure;
subplot (3,1,1);
plot (decimate(all_data (:,2),dr),decimate(all_data (:,3),dr));
xlabel('East Position(m)');
ylabel('Down Position(m)');
hold on;
quiver (decimate(all_data (:,2),dr),decimate(all_data (:,3),dr),decimate(all_data (:,4),dr),decimate(all_data (:,5),dr));
legend('Roll at position');
xlabel('East Position(m)');
ylabel('Down Position(m)');

%quiver plot pitch
subplot (3,1,2);
plot(decimate(all_data (:,2),dr),decimate(all_data (:,3),dr));
xlabel('North Position(m)');
ylabel('Down Position(m)');
hold on;
quiver (decimate(all_data (:,2),dr),decimate(all_data (:,3),dr),decimate(all_data (:,6),dr),decimate(all_data (:,7),dr));
legend('Pitch at position');
xlabel('North Position(m)');
ylabel('Down Position(m)');

hold off;
subplot (3,1,3);
plot (decimate(all_data (:,2),dr),decimate(all_data (:,1),dr));
xlabel('East Position(m)');
ylabel('North Position(m)');
% Plot Yaw on comet plot
hold on;
quiver (decimate(all_data (:,2),dr),decimate(all_data (:,1),dr),decimate(all_data (:,8),dr),decimate(all_data (:,9),dr));
legend('Yaw at position');
xlabel('East Position(m)');
ylabel('North Position(m)');
hold off;

% plot(Time_window , Abs_pos, Time_window , Abs_vel, Time_window , Abs_accel)
% legend('Position', 'Velocity' , 'Acceleration');
% xlabel('Time(s)');
% ylabel('m,m/s,m/s2 ');


figure;
subplot (2,2,1)
plot(Time_window(1:end-2,1), Abs_pos(1:end-2,1),Time_window(1:end-2,1), Abs_vel(1:end-1,1),Time_window(1:end-2,1), Abs_accel)
legend('Position', 'Velocity' , 'Acceleration');
xlabel('Time/s');
ylabel('m,m/s,m/s2 ');

subplot (2,2,2)
plot( Time_window(1:end-2,1), Abs_accel);
legend( 'Acceleration magnitude');
xlabel('Time/s');
ylabel('m/s2 ');

subplot (2,2,3)
plot(decimate(all_data (:,2),dr),decimate(all_data (:,1),dr));
xlabel('East Position(m)');
ylabel('North Position(m)');



subplot(2,2,4)
plot (Angles_time(:,7), Angles_time(:,4),Angles_time(:,7), Angles_time(:,5),Angles_time(:,7), Angles_time(:,6));
xlabel('Time (s)');
ylabel('Roll, Pitch, Yaw/rad');
legend('Roll', 'Pitch', 'Yaw');
% 
%figure;
% plot3(all_data(:,2),all_data(:,1), all_data(:,3), position_east_gps, position_north_gps, data_alt_gps);
% legend('Fusion Position', 'GPS');
% xlabel('East (m)');
% ylabel('North (m)');
% zlabel('Altitude (m)');
% axis equal;


quiv_ds_rate = 50; % downsample rate

hold off;


figure;
quiver (decimate(NEDT(1:end-1,2), quiv_ds_rate), ...
    decimate(NEDT(1:end-1,1), quiv_ds_rate), ...
    decimate(diff(NEDT(:,2)), quiv_ds_rate), ...
    decimate(diff(NEDT(:,1)), quiv_ds_rate));
legend('Velocity at position');
xlabel('East Position(m)');
ylabel('North Position(m)');

figure;

plot3(NEDT(1:end-2,2),(NEDT(1:end-2,1)),NEDT(1:end-2,3));
hold on;
quiver (decimate(NEDT(1:end-2,2), quiv_ds_rate), ...
    decimate(NEDT(1:end-2,1), quiv_ds_rate), ...
    decimate(50*accelE(1:end,:), quiv_ds_rate), ...
    decimate(50*accelN(1:end,:), quiv_ds_rate));
legend('Acceleration at position');
xlabel('East Position(m)');
ylabel('North Position(m)');
hold off;

figure;
accelE = movmean(diff(diff(NEDT(:,2))),100);
accelN = movmean(diff(diff(NEDT(:,1))),100);
% smooth_accels = movmean(all_data_and_accels,50);
plot(all_data(:,2), all_data(:,1));
hold on;
quiver (NEDT(1:end-2,2),NEDT(1:end-2,1), accelE(1:end,:), accelN(1:end,:)); 

legend('Fusion Position','accel direction');
xlabel('East (m)');
ylabel('North (m)');

axis equal;
hold off;


figure;
plot3(kal_x_stor(2, :), kal_x_stor(1, :), -kal_x_stor(3, :));
legend('3D osition');
xlabel('East Position(m)');
ylabel('North Position(m)');
zlabel('Height(m)');
