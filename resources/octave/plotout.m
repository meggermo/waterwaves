function plotout

load dim.out;
load crd.out; 
x = crd(:,1); z = crd(:,2); 
u = crd(:,5); v = crd(:,6);
N = sum(dim);
I = 1:N;
T = floor(size(crd,1)/N);
cla;
hold off;
for i=1:T,J=I+N*(i-1);
  cla;
  hold on;
  quiver(x(J), z(J), u(J), v(J));
  plot(x(J), z(J), 'ro-');  
  axis([0 45 -4 2]);
  pause(0.02);
end