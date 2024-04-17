clear
r = 0.25;
goal = [5 -10 griewank([-5 -10])];
start = [0 0 1.1];
path = OptimalPathGradient(start,goal);

numpts = 25;
step = floor(length(path)/numpts);
waypoints(1,:) = start;
waypoints(numpts+1,:) = path(length(path),:);
for i = 1:numpts
    waypoints(i,:) = path(i*step,:);
end

for i = 1:numpts
    [y(:,(i-1)*1000+1:i*1000), tau((i-1)*1000+1:i*1000)] = TPBVP2(waypoints(i,:),waypoints(i+1,:));
    tau((i-1)*1000+1:i*1000)= ((i-1))+tau((i-1)*1000+1:i*1000);
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