function [dy] = posDotInterceptor(t,y)
%this is the second (third) trial algo for navigation. It consists on this:
%given the fact that I have a certain velocity, try to intercept the target
%on its courrent course, assuming he will move in a straight line

%said velocity is, FOR EXAMPLE (FIRST TRIAL) such that in a second I would
%cover 10-20% of the (crurrent, for now) distance to target. Of course, there is a cap in the
%speed i can reach (say 5 ms/s)

%pos IS A 3dim vecotr, NOT 6-DIM LIKE IN PNG (control is made not in
%acceleration, but in speed)



%this piece of code is added for compatibily. It's not the best thing in
%the world...
global tgtPos tgtTimes pos currentStep;


lPos=pos(currentStep,:);
lTgtPos=tgtPos(currentStep,:);
lTargetVel=(tgtPos(currentStep,:)-tgtPos(currentStep-1,:))/(tgtTimes(currentStep)-tgtTimes(currentStep-1));
%end of compatibility block
%l (L) stands for local, in order not to use the same name as global
%variables

VEL_CAP=7;
K=1; %0.4

dist=norm(lPos-lTgtPos);
vel=min(K*dist,7);
lPredictedPos=lTgtPos+lTargetVel*dist/vel;
dist=norm(lPos-lPredictedPos);
v=(lPredictedPos-lPos)/dist;
vel=vel*v;
dy=vel';
end

