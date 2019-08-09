function [tgt]=targetTrajectory(t)
%this function should provide the position of the target, given the current
%time

tgt=zeros(1,3);
%uniform speed
v=[1 2 3];
vel=8.5;
v=v./norm(v);
tgt=t*v*vel;
tgt=tgt+[10 0 0];

%eg: target far away doing a circle
%  w=0.3; R=15;
%  tgt=zeros(1,3);
%  tgt(1)=10+R*cos(w*t);
%  tgt(2)=R*sin(w*t);
%  tgt(3)=0;

%eg: tgt uniformly accelerating
% a=0.7;
% tgt=zeros(1,3);
% tgt(1)=10;
% tgt(2)=0;
% tgt(3)=a*t^2/2;

%tgt uniformly accelerating but with capped speed
% a=0.7;
% tgt=zeros(1,3);
% tgt(1)=10;
% tgt(2)=0;
% vcap=2;
% tcap=vcap/a;
% if(t<tcap)
%     tgt(3)=a*t^2/2;
% else
%     tgt(3)=a*tcap^2/2+vcap*(t-tcap);
% end

%spiral
%  w=5*2*pi/15; R=3; vAx=2;
%  ax=[-3 2 6]; ax=ax/norm(ax);
%  n1=[2.1633    1.2245    0.6735];
%  n2=[0.8571   -2.1429    1.1429];
%  tgt=ax*vAx*t + R*(n1*cos(w*t)+n2*sin(w*t));
%  tgt=tgt./10;
%  tgt=tgt+[10 0 0];

end