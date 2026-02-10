u = all_data_vel (:,1:2);
v = all_data_vel(:,8:9);
z = 0;
    %join all data velocity and roll pitch yaw and equalise matrix size
for z = 1 :length(all_data_vel)    
    
    magu (z,:) = [norm(u(z,1:2))];

end
z = 0;
for z = 1 :length(all_data_vel)    
   
    magv (z,:) = [norm(v(z,1:2))];
end
figure;
CosTheta  = (dot(u,v,2))./(magu.* magv);
ThetaInDegrees = acosd(CosTheta);
plot (all_data_vel(:,10), ThetaInDegrees);
title('Angle between velocity and Yaw rotation');

xlabel('Time/s');
ylabel('Angle of attack/ degrees ');
z = 0;

figure;
z = 0;
dr = 50;

%comet
plot (all_data(:,2),all_data(:,1));
axis([-50 50 -50 50])
hold on;
quiver (decimate(all_data(:,2),dr), decimate(all_data(:,1),dr), decimate(all_data (:,8),dr), decimate( all_data (:,9),dr));
% legend('Yaw at position');

quiver (decimate(all_data(:,2),dr), decimate(all_data(:,1),dr), decimate(all_data_vel (:,2),dr), decimate( all_data_vel (:,1),dr));
hold off;
legend('Yaw at position','Velocity at position');

xlabel('East Position(m)');
ylabel('North Position(m)');

