clear
r = 0.25;
goal = [8.2 4.4 200*griewank([8.2 4.4])+1.1*r];
start = [0 0 1.1];
path = OptimalPathGradient(start,goal);

numpts = 100;
step = floor(length(path)/numpts);
waypoints(1,:) = start;
waypoints(numpts+1,:) = goal;
for i = 2:numpts
    
    waypoints(i,:) = path(i*step,:);
end