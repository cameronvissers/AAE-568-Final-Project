clear
r = 0.25;
% 10 6
% 4 -10
% -8.2 -4.4
% -6.6 9
% -5 -10
% 0 10
goal = [-8.2 -4.4 griewank([-8.2 -4.4])];
start = [0 0 1.1];
path = OptimalPathGradient(start,goal);
v = 15;

numpts = 40;
waypoints = zeros(numpts+1,12);
movedist = zeros(numpts,1);
movetime = zeros(numpts,1);
step = floor(length(path)/numpts);
waypoints(1,1:3) = start;
waypoints(numpts+1,1:3) = path(length(path),:);
for i = 2:numpts
     waypoints(i,1:3) = path(i*step,:);
     movedist(i-1) = sqrt((waypoints(i-1,1)-waypoints(i,1))^2+(waypoints(i-1,2)-waypoints(i,2))^2+(waypoints(i-1,3)-waypoints(i,3))^2);
     movetime(i-1) = movedist(i-1)/v;
     waypoints(i,4) = v*(waypoints(i,1)-waypoints(i-1,1))/movedist(i-1);
     waypoints(i,5) = v*(waypoints(i,2)-waypoints(i-1,2))/movedist(i-1);
     waypoints(i,6) = v*(waypoints(i,3)-waypoints(i-1,3))/movedist(i-1);
     if movedist(i-1) <0.05
         break
     end
end
%waypoints(numpts+1,4:6) = waypoints(numpts,4:6);
k = i-1;
for i = 1:k
    if i < 5 
        dt = 4;
    else
        dt = 1;
    end
    [y(:,(i-1)*1000+1:i*1000), tau((i-1)*1000+1:i*1000)] = TPBVP2(waypoints(i,:),waypoints(i+1,:),dt);
    if i >1
        tau((i-1)*1000+1:i*1000)= tau((i-1)*1000)+tau((i-1)*1000+1:i*1000);
    end
end
% g=9.81;
m=1;
% %I = [Ixx Iyy Izz]
%I = [0.0012 0.0012 0.002];
I = [1 1 1];
for i =1:length(y)
    vel(i) = sqrt(y(2,i)^2+y(4,i)^2+y(6,i)^2);
    u1(i) = min(300,max(-300,-y(18,i)/m));
    u2(i) = min(300,max(-300,-y(20,i)/I(1)));
    u3(i) = min(300,max(-300,-y(22,i)/I(2)));
    u4(i) = min(300,max(-300-y(24,i)/I(3)));
end

figure(2)
hold on
plot(waypoints(:,1),waypoints(:,2))

figure(1)
clf
plot3(waypoints(:,1),waypoints(:,2),waypoints(:,3))
hold on
plot3(y(1,:),y(3,:),y(5,:))

figure(3)
clf
plot(tau,-y(18,:),tau,-y(20,:),tau,-y(22,:),tau,-y(24,:))

figure(5)
clf
plot(tau,y(1,:),tau,y(3,:),tau,y(5,:))
title('Position over time')
xlabel('t seconds')
ylabel('Cartesian coordinate (m)')


figure(7)
clf
plot(tau,vel)
title('Quadcopter Total velocity over time')
xlabel('t seconds')
ylabel('Velocity m/s')

figure(8)
clf
sgtitle('Control Inputs')
subplot(2,2,1)
plot(tau,u1)
xlabel('t seconds')
ylabel('u1 input to $\ddot{z}$','interpreter','Latex')

subplot(2,2,2)
plot(tau,u2)
xlabel('t seconds')
ylabel('u2 input to $\ddot{\phi}$','interpreter','Latex')

subplot(2,2,3)
plot(tau,u3)
xlabel('t seconds')
ylabel('u3 input to $\ddot{\theta}$','interpreter','Latex')

subplot(2,2,4)
plot(tau,u4)
xlabel('t seconds')
ylabel('u4 inputs to $\ddot{\psi}$','interpreter','Latex')

