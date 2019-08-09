clear all;
%IMPORTANT CHOICE:I want every single position and direction to be a 1x3 ROW vector
%this must then be transposed inside the differential equation solver!
algoType=2;


startTime=0;
endTime=14;
timeStep=0.001;
nStep=(endTime-startTime)/timeStep; %start time is always 0

global tgtPos tgtTimes pos currentStep;
if algoType == 1
    global accNorm; %pre-allocate the variables, otherwise the solver will be much slower
end
tgtPos=ones(nStep,3);
tgtTimes=ones(nStep,1);
pos=ones(nStep,3);
if algoType == 1
    accNorm=ones(nStep,1);
end

currentStep=1;

%this is the first handmade integration step, in order not to have
%problems. See the function posDot.
switch algoType
    case 1
        y0=[0;0;0;1;0;1]; %start in the origin moving along x-z byseptric(?)
    case 2
        y0=[0;0;0];
end
    
aux=targetTrajectory(startTime);
tgtPos(currentStep,:)=aux;
tgtTimes(currentStep)=startTime;
pos(currentStep,:)=y0(1:3);
if algoType == 1 
    accNorm(currentStep)=0;
end
currentStep=currentStep+1;

switch algoType
    case 1
        y0(1:3)=y0(1:3)+y0(4:6)*1e-5;
        startTime=startTime+1e-5;
        currenty=y0;
    case 2
        currenty=y0;
end
        


assert(currentStep==2);
for currentStep=2:nStep
    switch algoType
        case 1
            dy=posDotPngTrue3D(currentStep*timeStep,currenty); %for first png algorithm
        case 2
            
            dy=posDotInterceptor(currentStep*timeStep,currenty);
    end
    currenty=currenty+dy*timeStep;
end


%options= odeset('RelTol',1e-5,'AbsTol',1e-5);
%[t, y]=ode23(@posDot,[startTime 7],y0,options); DONT USE THE NATIVE SOLVER



