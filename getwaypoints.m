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

figure(1)
plot3(waypoints(:,1),waypoints(:,2),waypoints(:,3))