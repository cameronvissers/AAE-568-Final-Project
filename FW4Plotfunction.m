xmin = -10;
xmax = 10;

f_size=14; % font size

for i = 1:101;
x1(i) = (xmax-xmin)/100*(i-1)+xmin;
x2(i) = (xmax-xmin)/100*(i-1)+xmin;
end

for i = 1:101
    for j =  1:101
        f(i,j) = 200*griewank([x1(j) x2(i)]);
    end
end

[X, Y] = meshgrid(x1,x2);

a = figure(4);
a = surf(X,Y,f);
set(gca,'Fontsize',f_size)
xlabel('x','Fontsize',f_size);
ylabel('y','Fontsize',f_size);
zlabel('h','Fontsize',f_size);
box on
hold on

[Xg, Yg] = meshgrid(x1/.0025,x2/.0025);

b = figure(6);
b = surf(Xg,Yg,f);
set(gca,'Fontsize',f_size)
xlabel('x','Fontsize',f_size);
ylabel('y','Fontsize',f_size);
zlabel('h','Fontsize',f_size);
box on
hold on

c = figure(2);
c = contour(X,Y,f,20);
colorbar
set(gca,'Fontsize',f_size)
xlabel('x','Fontsize',f_size);
ylabel('y','Fontsize',f_size);
zlabel('h','Fontsize',f_size);
box on;
hold on
