function [dy] =posDotPngTrue3D(t,y)
    global tgtPos tgtTimes pos currentStep accNorm;
    aux=targetTrajectory(t);
    tgtPos(currentStep,:)=aux;
    tgtTimes(currentStep)=t;
    pos(currentStep,:)=(y(1:3))'; %BEWARE:Y COMES AS A COLOUMN VECTOR(?)
    
    %WE MUST BEGIN by having ALREADY taken a SMALL integration step before
    %beginning; Initial conditions must be given not from time 0 but from, eg,
    %time 1e-5. In this step, we suppose that velocity is constant.
    currTgtPos=tgtPos(currentStep,:);
    lastTgtPos=tgtPos(currentStep-1,:);
    
    currMyPos=pos(currentStep,:);
    lastMyPos=pos(currentStep-1,:);
    
    currTime=tgtTimes(currentStep);
    lastTime=tgtTimes(currentStep-1);
    [lambdaDot, v, closingVel,LOS]=LOSVariations([lastTgtPos;currTgtPos],[lastTime currTime],[lastMyPos; currMyPos],[],2);
    
    if(~isnan(v))
        normalAcc=300*lambdaDot*closingVel*v;
    else
        normalAcc=[0 0 0];
    end
    forwardAcc=tanAccControl(closingVel,LOS);
    acc=normalAcc+forwardAcc;
    accNorm(currentStep)=norm(acc);
    dy=zeros(6,1);
    dy(1)=y(4); %x dot=Vx;
    dy(2)=y(5);
    dy(3)=y(6);
    
    dy(4)=acc(1);
    dy(5)=acc(2);
    dy(6)=acc(3);
    
    %assert(acc(3)==0); 
    
%     closingVel
%     currentStep
%     disp("-----------------");

     %input("--------------------------------");
    currentStep=currentStep+1;
end