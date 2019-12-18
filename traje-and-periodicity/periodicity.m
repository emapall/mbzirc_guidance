clear all;
close all;
global x;
global y;
global t;
global l;
global ERR_TOL_COEFF;
global HALF_TIME_TOL;
HALF_TIME_TOL=3; %the tolerance with wich we consider that two points are "near"
ERR_TOL_COEFF=2;
a=10;
t=0:0.1:15*2;
%starts from upper right corner, spins cw, 
%10 secs==1 spin +right arm and low r- up l arm
%15 secs==more than 2 complete rounds
samplesPerWin=25; 


l=length(t);
devi=a/20; 
%x is comprised between -a and a, y between a/sqrt(2) and -a/sqrt(2)
[x,y]=traj(a,t);
noisex=normrnd(0,devi,1,l);
noisey=normrnd(0,devi,1,l);
x=x+noisex;
y=y+noisey;

% hold on
% for i=1:l
%     plot(x(i),y(i),"rx");
%     pause(0.1);
% end

[err,p]=prelimAnalysis(samplesPerWin); 
plot(err,"rx");
hold on
for i=p
    plot(i,err(i),"ko");
end

figure;
subplot(2,1,2);
plot(x,y,"k+");
global xl;
global yl;
yl=ylim;
xl=xlim;
global sp;
sp=subplot(2,1,1);
ylim(yl);
xlim(xl);
xlim manual;
ylim manual;

dists=zeros(1,2,1);
cci=1; %comulative dandiate index, for readability
iGiro=0; %il primo giro è lo 0, dento il while si parte da 1

while(cci<=length(p))
    candidatiGiro=find(abs(p-p(cci))<=2*HALF_TIME_TOL);
    cci=cci+length(candidatiGiro);
    iGiro=iGiro+1;
    
    for i=1:length(candidatiGiro) %for candidato=p(candidatiGiro)
       [temp, idxCorrels]=mobileWindowDist(samplesPerWin,iGiro,p(candidatiGiro(i)));
       dists(i,length(temp),iGiro)=-422;
       length(temp)
       dists(i,1:length(temp),iGiro)=temp;
    end
end

figure;


function [x,y]=traj(a,t)
    x=a*2*sin(t)./(sin(t).^2+1.);
    y=x.*cos(t);
end


function [err,periodIdxCandidates]=prelimAnalysis(nWin)
%calcola la distanza (?) tra 2 finestre mobili
%per ora è la SOMMA delle distanze euclidee tra i punti
%err(1) è sempre 0
    global x;
    global y;
    global l;
    global ERR_TOL_COEFF;
    w0x=x(1:nWin);
    w0y=y(1:nWin);
    minPeriod=-1;
    err=-10*ones(1,l);
    for i=2:l-nWin
        winx=x(i:i+nWin-1);
        winy=y(i:i+nWin-1);
        err(i)=distanceBetweenTwoWins((w0x-winx),(w0y-winy));
    end
    periodIdxCandidates=find(err<=ERR_TOL_COEFF*min(err(2:l-nWin)) & err>=0);  
    periodIdxCandidates=periodIdxCandidates(2:length(periodIdxCandidates));
    %LEAVE the first one apart;
    while(periodIdxCandidates(1)<=ERR_TOL_COEFF*min(err(2:l-nWin)))
        periodIdxCandidates=periodIdxCandidates(2:end);
    end
    
    %TODO:controllare che il minimo non sia in posizione 2, quando
    %l'errore deve ancora salire
end

function [err,idxCorrels]=mobileWindowDist(nWin,nGiro,idxStart)
    %resittuisce un vettore di distanza tra le finestre mobili 
    %"corrispondenti".

    %il problema grosso è trovare cosa siano queste "corrispondenze": come
    %criterio, scielgo il seguente (discutibile): per ogni punto, il punto di
    %traiettoria che gli dovrebbe essere "corrispondente" al giro prima è
    %quello che più si avvicina col suo timestamp temporale al timestamp del
    %giro dopo, meno il timestamp del periodo

    global x;
    global y;
    global t;
    global l;
    

    %T=t(idxStart)/nGiro;
    T=t(idxStart); 
    %per sciegliere tra le 2, bisogna vedere se vogliamo fare il confronto
    %tra il PRIMO giro e quiesto giro, o tra il giro prima e questo giro
    %(cfr findNearestTemporalCorrespondant)
    
    %we cannot know in advance how many samples we will have for THIS turn
    err=-42*ones(1,ceil(idxStart/nGiro) );

    idx=idxStart;
    idxCorr=-1;
    idxCorrels=err;

    while(t(idx)<T*(nGiro+1)/nGiro & idx+nWin-1<=l)
        idxCorr=findNearestTemporalCorrespondant(idx,T);
        xCur=x(idx:idx+nWin-1);
        xPrev=x(idxCorr:idxCorr+nWin-1);
        yCur=y(idx:idx+nWin-1);
        yPrev=y(idxCorr:idxCorr+nWin-1);
        err(idx-idxStart+1)=distanceBetweenTwoWins((xCur-xPrev),(xCur-xPrev));
        
        subplotWindows(x(idx:idx+nWin-1),y(idx:idx+nWin-1),x(idxCorr:idxCorr+nWin-1),y(idxCorr:idxCorr+nWin-1));
        %pause(0.05);
        
        idxCorrels(idx)=idxCorr;
        idx=idx+1;
    end
end


function [idxPrev]=findNearestTemporalCorrespondant(idxCurr,T)
    global t;
    a=t(1:idxCurr); %avoid successive loops
    a=abs(a- (a(idxCurr)-T) );
    idxPrev=find(a==min(a));
    idxPrev=idxPrev(1);
end
   
function [d]=distanceBetweenTwoWins(a,b)
%return the distance between two vectors, according to the most suitable
%method we can find
    d=sum(sqrt( a.^2+b.^2 ));
end

function[]=subplotWindows(xcur,ycur,xpast,ypast)

    global xl;
    global yl;

    plot(0,0);
    hold on;
    ylim(yl);
    xlim(xl);
    xlim manual;
    ylim manual;
    plot(xcur,ycur,"rx");
    plot(xpast,ypast,"ks");
    hold off;
end