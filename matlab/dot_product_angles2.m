u = 100*(diff(all_data (:,1:2)));
v = all_data(1:end-1,8:9);
A = zeros((length(UT)-1), 1);
uc = [u A];
vc = [v A];
z = 0;
    %join all data velocity and roll pitch yaw and equalise matrix size
    magu = zeros;
    magv = zeros;
 
magu = sqrt((u(:,1).^2) + (u(:,2).^2));
magv = sqrt((v(:,1).^2) + (v(:,2).^2));  
    
% figure;
% cosTheta = []
% ThetaInDegress = []
% for b = 1 :length(uc)
%     
% CosTheta(b)  = (dot(u(b,:),v(b,:),2))./(magu(b).* magv(b));

plot3 (alldata(:,3) , alldata(:,2), alldata(:,1)) 
% ThetaInDegreesA(b) = acosd(CosTheta(b));
% 
% end

CP = (cross(uc,vc,2));
CPmag = sqrt(CP(:,3).^2);
ThetaInDegrees = [];
for a = 1 :length(uc)

ThetaInDegrees(a) = atan2d(CPmag(a,1),dot(uc(a,1),vc(a,1),2));
end
figure
plot (all_data(1:end-1,10), movmean(ThetaInDegrees,50));
title('Angle between velocity and Yaw rotation');
xlabel('Time/s');
ylabel('Angle of attack/ degrees ');

% figure
% 
% plot (all_data(1:end-1,10), ThetaInDegreesA);
% title('Angle between velocity and Yaw rotationA');
% xlabel('Time/s');
% ylabel('Angle of attack/ degrees ');
% z = 0;

% figure;
% z = 0;
% dr = 1;
% comet (all_data(:,2),all_data(:,1));
% % axis([-5 15 -5 10])
% hold on;
% quiver (decimate(all_data(1:end-1,2),dr), decimate(all_data(1:end-1,1),dr), decimate(v (:,1),dr), decimate( v(:,2),dr));
% % legend('Yaw at position');
% 
% quiver (decimate(all_data(1:end-1,2),dr), decimate(all_data(1:end-1,1),dr), decimate(u(:,2),dr), decimate(u(:,1),dr));
% hold off;
% legend('Yaw at position','Velocity at position');
% 
% xlabel('East Position(m)');
% ylabel('North Position(m)');
% 
figure;
plot (all_data(:,2),all_data(:,1));
hold on;
quiver (decimate(all_data(1:end-1,2),dr), decimate(all_data(1:end-1,1),dr), decimate(u(:,2),dr), decimate(u(:,1),dr));
hold off;
% 
% figure;
% plot (all_data(:,2),all_data(:,1));
% hold on;
% quiver (decimate(all_data(1:end-1,2),dr), decimate(all_data(1:end-1,1),dr), decimate(v(:,1),dr), decimate(v(:,2),dr));
% hold off;
