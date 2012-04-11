function sim_select_test(n)

if nargin<1, n=1000; end
x=10*rand(n,2)-5;
[y,z]=sim_select(x,x,[10 10]);
close(gcf);
plot(x(:,1),x(:,2),'.',y(:,1),y(:,2),'.',z(:,1),z(:,2),'o');
grid