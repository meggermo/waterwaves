function waveplot

load dim.out;
load crd.out; x = crd(:,1); z = crd(:,2);
N=sum(dim);
I=1:N;
T=floor(size(crd,1)/N);
for i=T-10:T-2,J=I+N*(i-1);
  plot(x(J),z(J),'o-');
  axis([0 40 -4 1]);
  pause(1);
end