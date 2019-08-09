function [a] = tanAccControl(closingVel,LOS)
%a simple pid for acceleration along line of sight
%acceleration is capped to +- 3
kVel=0.3; %1 of acceleration for every 3 m/s of error upon the speed
v0=1.5; 
kDist=0.1;
d0=10;
a=0;
dist=norm(LOS);
n=LOS/norm(LOS);
assert(norm(LOS)~=0); %otherwise....

a=kVel*(v0-closingVel);

if(dist>=d0)
    a=a+kDist*dist;
end

a=a*n;

end

