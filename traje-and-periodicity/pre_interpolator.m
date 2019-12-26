clear all;
close all;
global x;
global y;
global t;
global l;
global ERR_TOL_COEFF;
global HALF_TIME_TOL;
global TOTAL_SIM_TIME;
global MIN_TIME_INC;
global MAX_TIME_INC;
global SIM_TIME_STEP;
HALF_TIME_TOL=1; %the tolerance with wich we consider that two points are "near"
ERR_TOL_COEFF=2;
TOTAL_SIM_TIME=31.4;
SIM_TIME_STEP=0.1;
MIN_TIME_INC=SIM_TIME_STEP/2;
MAX_TIME_INC=SIM_TIME_STEP*1.5;
TIME_VARIANCE_COEFF=SIM_TIME_STEP/0.2;
a=100;
t=0:SIM_TIME_STEP:TOTAL_SIM_TIME;
%starts from upper right corner, spins cw, 
%10 secs==1 spin +right arm and low r- up l arm
%15 secs==more than 2 complete rounds
samplesPerWin=8; 


l=length(t);
devi=a/7; 
%x is comprised between -a and a, y between a/sqrt(2) and -a/sqrt(2)
[x,y]=traj(a,t);
noisex=normrnd(0,devi,1,l);
noisey=normrnd(0,devi,1,l);
x=x+noisex;
y=y+noisey;
%% 
[err,p]=prelimAnalysis(samplesPerWin); 
%% GRAFICI
subplot(2,2,1);
hold on;
subplot(2,2,3);
hold on;
for i=1:l+1-samplesPerWin
    subplot(2,2,1);
    plot(t(i),err(i),"r.");
    %TODO FIX: stampare coi tempi, non con gli indici
    subplot(2,2,3);
    plot(x(i),y(i),"ko");
    %pause(0.03);
end

for i=l-samplesPerWin+2:l
    plot(x(i),y(i),"o");
end
subplot(2,2,1);
hold on;
for i=p
    plot(t(i),err(i),"ys");
end

%% ANALISI GIRI

global tApproxGiri;
tApproxGiri=zeros(1,2);
tempidx=1;

while(tempidx<=length(p))
    auxvar=find(abs(t(p)-t(p(tempidx)))<=2*HALF_TIME_TOL);
    %TODO: FIX: il criterio dev'essere temporale, non sugli indici
    tApproxGiri(end+1)=sum(t(p(auxvar)))/length(auxvar);
    tempidx=tempidx+length(auxvar);
end
tApproxGiri=tApproxGiri(2:end); %uno zero lo lascio che torna comodo fra poche righe
nGiri=length(tApproxGiri)-1; %NUMERO DI GIRI COMPLETATI 
%ha due zeri all'inizio per inizializzazione 
tApproxMid=tApproxGiri;
tApproxMid(2:end)=(tApproxGiri(2:end)+tApproxGiri(1:end-1))/2;
%%
pause(1.5);

%% ANALISI SPAZIALE

idxTest=l;
distVicini=sqrt((x-x(idxTest)).^2 +((y-y(idxTest)).^2));
distVicini(idxTest)=max(distVicini);
minDist=min(distVicini);

tolprova5=3;
subplot(2,2,2);
hold on;
for(toll_coeff=1:0.05:7)
    %[toll_coeff sum(distVicini<=minDist*toll_coeff)]
    subplot(2,2,2);
    plot(toll_coeff,sum(distVicini<=minDist*toll_coeff),"r.");
    
    subplot(2,2,3);
    if(toll_coeff==round(toll_coeff))
        circle(x(idxTest),y(idxTest),toll_coeff*minDist);
    end
    if(sum(distVicini<=minDist*toll_coeff)<=5 && sum(distVicini<=minDist*toll_coeff)>=3)
        tolprova5=toll_coeff;
    end
end

tolprova=3;

pastCandSpaz=find(distVicini<=minDist*tolprova);
subplot(2,2,1);
hold on;
for i=pastCandSpaz
    if(i<=l-samplesPerWin+1)
        plot(t(i),err(i),"bx");
    else
        plot(t(i),min(err(2:end))*2.5,"gs");
    end
end

giroPastCand=-1+zeros(1,length(pastCandSpaz));
for i=1:length(pastCandSpaz)
    for j=1:nGiri
        if(tApproxGiri(j) < t(pastCandSpaz(i)) && t(pastCandSpaz(i))<=tApproxGiri(j+1))
            giroPastCand(i)=j-1;
        end
    end
end

%%
%POSSIAMO SCIEGLIERE IL CRITERIO TEMPORALE COME BUON CRITERIO! QUELLO
%SPAZIALE, A MENO DI NON PRENDERE GRANDI RAGGI, NON FUNZIONA BENISSIMO
%% ANALISI TEMPORALE: SE  varianza sul periodo è troppo alta, allora non si usa
%if (varianza(tApproxGiri)<=soglia)
T=mean(diff(tApproxGiri));
%percGiro=(t(idxTest)-tApproxGiri(end))/T;
%ti piacissi: è un criterio non buono; Usiamo il seguente:
%trova il corrispondente spaziale più vicino all'ultimo campione
%registrato, per ogni giro passato ; Per ognuno
% di quei campioni corrispondenti calcola la percentuale di
%completamento giro, e ci fa la media sopra
percGiro=-1+zeros(1,nGiri);
pastCandQuasiDist=-1*ones(1,nGiri);

subplot(2,2,3);
for iGiro=1:nGiri
    idxPrev=findNearestSpatialCorrespondant(x(idxTest),y(idxTest),iGiro);
    pastCandQuasiDist(iGiro)=idxPrev;    
    plot(x(idxPrev),y(idxPrev),"rx");
    percGiro(iGiro)=(t(idxPrev)-tApproxGiri(iGiro))/(tApproxGiri(iGiro+1)-tApproxGiri(iGiro));
    assert(percGiro(iGiro)>=0 &&percGiro(iGiro)<=1);   
end
   
percLast=mean(percGiro);
pastCandTemp=zeros(1,nGiri); %se un giro è appena appena completato

for i=1:nGiri%-1
    tCorrispIdeale=tApproxGiri(i)*(1-percLast)+tApproxGiri(i+1);
    auxvar=abs(t-tCorrispIdeale);
    pastCandTemp(i)= min(find(auxvar==min(auxvar)));
    %va gestito il caso in cuici siano 2 elementi che siano ESATTAMENTE il min
    %può capitare sia quasi 0, occhio agli errori numerici!
    subplot(2,2,1);
    if(pastCandTemp(i)<=l-samplesPerWin-1)
        plot(pastCandTemp(i),err(pastCandTemp(i)),"kh");
    end
    subplot(2,2,3);
    plot(x(i),y(i),"rh");
    
end
%% costruzione insieme di punti
pastCand=pastCandTemp;
xWin=x(end-samplesPerWin+1:end);
yWin=y(end-samplesPerWin+1:end);

subplot(2,2,4);
hold on;
for i=pastCand
    auxvar=x(i-samplesPerWin+1:i+samplesPerWin-1);
    xWin(end+1:end+length(auxvar))=auxvar;
    auxvar=y(i-samplesPerWin+1:i+samplesPerWin-1);
    yWin(end+1:end+length(auxvar))=auxvar;
    
    plot(x(i-samplesPerWin+1:i+samplesPerWin-1),y(i-samplesPerWin+1:i+samplesPerWin-1),"x:");
    pause(1);
end

%% misto vs temp
g=figure;
hold on;
plot(x,y,"k.");
plot(x(idxTest),y(idxTest),"bo");
plot(x(pastCandTemp),y(pastCandTemp),"rx:");
plot(x(pastCandQuasiDist),y(pastCandQuasiDist),"gs:");

for i=pastCandTemp
    auxvar=x(i-samplesPerWin+1:i+samplesPerWin-1);
    xWin(end+1:end+length(auxvar))=auxvar;
    auxvar=y(i-samplesPerWin+1:i+samplesPerWin-1);
    yWin(end+1:end+length(auxvar))=auxvar;
    
    plot(x(i-samplesPerWin+1:i+samplesPerWin-1),y(i-samplesPerWin+1:i+samplesPerWin-1),"r:");
    pause(1);
end

for i=pastCandQuasiDist
    auxvar=x(i-samplesPerWin+1:i+samplesPerWin-1);
    xWin(end+1:end+length(auxvar))=auxvar;
    auxvar=y(i-samplesPerWin+1:i+samplesPerWin-1);
    yWin(end+1:end+length(auxvar))=auxvar;
    
    plot(x(i-samplesPerWin+1:i+samplesPerWin-1),y(i-samplesPerWin+1:i+samplesPerWin-1),"g:");
    pause(1);
end


%% funzioni
function circle(cx,cy,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.005:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(cx+xp,cy+yp);
end

function [x,y]=traj(a,t)
    %lookup table per il tempo: il tempo giro è variabile
    global MAX_TIME_INC;
    global MIN_TIME_INC;
    figure;
    plot(t);
    hold on;
    for i=2:length(t)
        t(i)=rand*(MAX_TIME_INC-MIN_TIME_INC)+MIN_TIME_INC+t(i-1);
    end
    plot(t);
    hold off;
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
    err=-10*ones(1,l-nWin+1);
    for i=2:l-nWin+1
        winx=x(i:i+nWin-1);
        winy=y(i:i+nWin-1);
        err(i)=distanceBetweenTwoWins((w0x-winx),(w0y-winy));
    end
    periodIdxCandidates=find(err<=ERR_TOL_COEFF*min(err(2:l-nWin)) & err>=0);  
    %LEAVE the first one apart;
    if(periodIdxCandidates(1)==2)
        while(err(periodIdxCandidates(1)+1)<=ERR_TOL_COEFF*min(err(2:l-nWin)))
            err(periodIdxCandidates(1))
            periodIdxCandidates=periodIdxCandidates(2:end);
        end
        periodIdxCandidates=periodIdxCandidates(2:end);
    end

    %TODO:controllare che il minimo non sia in posizione 2, quando
    %l'errore deve ancora salire
end

% function [idxPrev]=findNearestTemporalCorrespondant(idxCurr,T)
%     FUNZIONE DISMESSA
%     global t;
%     a=t(1:idxCurr); %avoid successive loops
%     a=abs(a- (a(idxCurr)-T) );
%     idxPrev=find(a==min(a));
%     idxPrev=idxPrev(1);
% end
%% funzioni 2
function [idxPrev]=findNearestSpatialCorrespondant(xP,yP,giro)
%trova il punto più vicono al punto in input
    global x;
    global t;
    global y;
    global tApproxGiri;
    xlocal=x;
    ylocal=y;
    inizio=findNearestTemporalTimestap(tApproxGiri(giro));
    fine=findNearestTemporalTimestap(tApproxGiri(giro+1));
    xlocal=x(inizio:fine);
    ylocal=y(inizio:fine);
    
    aux=sqrt((xlocal-xP).^2 +(ylocal-yP).^2);
    idxPrev=find(aux==min(aux));
    idxPrev=idxPrev(1)+inizio-1;
end

function [idx]=findNearestTemporalTimestap(tP)
%ritorna l'indice del campione il cui timestamp più si avvicina a quel
%tempo tP
    global t;
    idx=find(abs(t-tP)==min(abs(t-tP)));
    idx=idx(1); %non si sa mai
end

function [d]=distanceBetweenTwoWins(a,b)
%return the distance between two vectors, according to the most suitable
%method we can find
    d=sum(sqrt( a.^2+b.^2 ));
end