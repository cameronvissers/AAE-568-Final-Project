function path = OptimalPathGradient2(start,goal) 

FW4Plotfunction
r = 0.5;
%goal = [-8.2 -4.4 griewank([-8.2 -4.4])+1.1*r];
%start = [0 0 1.1];
%stepsize
alpha = 0.001;


xmin = -10;
xmax = 10;

for i = 1:101
    x1(i) = (xmax-xmin)/100*(i-1)+xmin;
    x2(i) = (xmax-xmin)/100*(i-1)+xmin;
end

%defines obstacles repulsion and gradient mesh
for i = 1:101
    for j =  1:101
        Krepulsive(i,j) = griewank([x1(j) x2(i)]);
        [dFdx(i,j), dFdy(i,j)] = gradgriewank(x1(j),x2(i));
    end
end

[X, Y] = meshgrid(x1,x2);

current = start;
%loop until you reach goal or for 10 steps
for k=1:30000
    %check if you reached the goal
    dg = gxyz(goal,current);
    if dg <=0.1
        break
    end
    
    %generate gradient steps from goal for x y and z
    for ii = 1:3
        sumdF(ii) = (current(ii)-goal(ii))/gxyz(current,goal);
    end
   
    %locate nearest obstacle
    dmin = 999;
    ob = [0 0 0];
    for i = 1:101
        for j =  1:101            
            di(i,j) = gxyz(current,[x1(j) x2(i) Krepulsive(i,j)]); 
            if di(i,j) <dmin
                dmin = di(i,j);
                ob = [x1(j) x2(i) Krepulsive(i,j)];
            end
        end
    end
    sumdF(1) = sumdF(1)+(ob(1)/2000+sin(ob(1))*cos(ob(2)/sqrt(2)))*(current(1)-ob(1))/dmin*ob(3);
    sumdF(2) = sumdF(2)+(ob(2)/2000+1/sqrt(2)*cos(ob(2))*sin(ob(2)/sqrt(2)))*(current(2)-ob(2))/dmin*ob(3);
    sumdF(3) = sumdF(3)+(current(3)-ob(3))*ob(3)/dmin;
    path(k,:) = current;
    for ii =1:3
        dstep(k,ii) = -alpha*sumdF(ii);
        
        current(ii) = current(ii)+dstep(k,ii);

    end
end
%project onto the surface
for i= 1:length(path)
    path(i,3) = max(path(i,3),200*griewank([path(i,1) path(i,2)]));
end
figure(4)
hold on
plot3(path(:,1),path(:,2),path(:,3),'-c',LineWidth=3)
end

function [dFdx dFdy] = gradgriewank(x,y)
    dFdx = 200*(x/2000+sin(x)*cos(y/sqrt(2)));
    dFdy = 200*(y/2000+cos(x)*sin(y/sqrt(2))/sqrt(2));
end

function dist = gxyz(goal,current)
   dist = sqrt((goal(1)-current(1))^2+(goal(2)-current(2))^2+(goal(3)-current(3))^2);
end

function dist = fxyz(goal,current)
   dist = sqrt((goal(1)-current(1))^2+(goal(2)-current(2))^2);
end