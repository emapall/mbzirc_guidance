clear all;
close all;
global x;
global y;
global t;
global l;
global ERR_TOL_COEFF;
global HALF_TIME_TOL;
HALF_TIME_TOL=3;
ERR_TOL_COEFF=2;
a=10;
t=1:0.1:15*2;
%starts from upper right corner, spins cw, 
%10 secs==1 spin +right arm and low r- up l arm
%15 secs==more than 2 complete rounds

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
[err,p]=decider(25); 
plot(err,"rx");

candidati=zeros(1);
iGiro=1; %il primo giro è lo 0

candidatiGiro=find(p<=p(1)+2*HALF_TIME_TOL);

candidati(1,length(candidatiGiro))=0;
candidati(1,:)=p(candidatiGiro);



function [x,y]=traj(a,t)
    x=a*2*sin(t)./(sin(t).^2+1.);
    y=x.*cos(t);
end


function [err,periodIdxCandidates]=decider(nWin)
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
    err=zeros(1,l);
    for i=2:l-nWin
        winx=x(i:i+nWin-1);
        winy=y(i:i+nWin-1);
        err(i)=sum(sqrt( (w0x-winx).^2+(w0y-winy).^2 ));
    end
    periodIdxCandidates=find(err<=ERR_TOL_COEFF*min(err(2:l-nWin)));  
    periodIdxCandidates=periodIdxCandidates(2:length(periodIdxCandidates));
    %LEAVE the first one apart;
    while(periodIdxCandidates(1)<=ERR_TOL_COEFF*min(err(2:l-nWin)))
        periodIdxCandidates=periodIdxCandidates(2:end);
    end
    
    %TODO:controllare che il minimo non sia in posizione 2, quando
    %l'errore deve ancora salire
end

function [err]=mobileWindowDist(nWin,nGiro,idxStart)
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

T=t(idxStart);
if(l>(nGiro+1)/nGiro*idxStart-1+nWin-1)
    