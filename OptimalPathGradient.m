FW4Plotfunction
r = 0.25;
goal = [8.2 4.4 200*griewank([8.2 4.4])+1.1*r];
start = [0 0 1.1];
%stepsize
alpha = 0.05;


xmin = -10;
xmax = 10;

for i = 1:101
    x1(i) = (xmax-xmin)/100*(i-1)+xmin;
    x2(i) = (xmax-xmin)/100*(i-1)+xmin;
end

%defines obstacles repulsion and gradient mesh
for i = 1:101
    for j =  1:101
        Krepulsive(i,j) = 200*griewank([x1(j) x2(i)]);
        [dFdx(i,j), dFdy(i,j)] = gradgriewank(x1(j),x2(i));
    end
end

[X, Y] = meshgrid(x1,x2);

current = start;
%loop until you reach goal or for 10 steps
for k=1:60000
    %check if you reached the goal
    dg = fxyz(goal,current);
    if dg <=0.1
        break
    end
    
    %generate gradient steps from goal for x y and z
    for ii = 1:3
        sumdF(ii) = (current(ii)-goal(ii))/fxyz(current,goal);
    end

    %check repulsive forces from each obstacle
    for i = 1:101
        for j =  1:101
            %check distance from each obstacle
            di(i,j) = fxyz(current,[x1(j) x2(i) Krepulsive(i,j)]); 
            
            if di(i,j) <= r
                Frepulsion(i,j) = Krepulsive(i,j)*(1/di(i,j)-1/r);
            elseif di(i,j) > 1.5*r
                Frepulsion(i,j) = 0;
            else
                Frepulsion(i,j) = Krepulsive(i,j)/di(i,j)^2;
            end

            sumdF(1) = sumdF(1)+(current(1)-x1(j))/di(i,j)*Frepulsion(i,j);
            sumdF(2) = sumdF(2)+(current(2)-x2(i))/di(i,j)*Frepulsion(i,j);
            sumdF(3) = sumdF(3)+(current(3)-Krepulsive(i,j))/di(i,j)*Frepulsion(i,j);
        end
    end
    path(k,:) = current;
    for ii =1:3
        dstep(k,ii) = -alpha*sumdF(ii);
        if k > 1
            if ii == 3
                lim = 20*alpha;
            else
                lim = alpha;
            end
            if abs(dstep(k,ii))>lim
                dstep(k,ii) = dstep(k-1,ii);
            end
        end
        
        current(ii) = current(ii)+dstep(k,ii);
        if ii==3
            current(ii) = max(current(ii),200*griewank([current(1) current(2)])+r);
        end
    end
end

figure(4)
hold on
plot3(path(:,1),path(:,2),path(:,3),'-c',LineWidth=3)

function [dFdx dFdy] = gradgriewank(x,y)
    dFdx = 200*(x/2000+sin(x)*cos(y/sqrt(2)));
    dFdy = 200*(y/2000+cos(x)*sin(y/sqrt(2))/sqrt(2));
end

function dist = fxyz(goal,current)
   dist = sqrt((goal(1)-current(1))^2+(goal(2)-current(2))^2+(goal(3)-current(3))^2);
end