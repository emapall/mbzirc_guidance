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
global samplesPerWin
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
fig0=figure(); %contiene grafico errore e numero di vicini in funz del raggio
subplot(2,1,1);
hold on;
plot(t(1:l+1-samplesPerWin),err(1:l+1-samplesPerWin),"r.");
for i=p
    plot(t(i),err(i),"ys");
end

fig1=figure(); %traiettoria con punti individuati dai vari criteri
hold on;
plot(x(1:l+1-samplesPerWin),y(1:l+1-samplesPerWin),"k.");
for i=l-samplesPerWin+2:l
    plot(x(i),y(i),"o");
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
idxTest=l;

% %% ANALISI SPAZIALE
% 
% 
% distVicini=sqrt((x-x(idxTest)).^2 +((y-y(idxTest)).^2));
% distVicini(idxTest)=max(distVicini);
% minDist=min(distVicini);
% 
% tolprova5=3;
% subplot(2,1,2);
% hold on;
% for(toll_coeff=1:0.05:7)
%     %[toll_coeff sum(distVicini<=minDist*toll_coeff)]
%     subplot(2,1,2);
%     plot(toll_coeff,sum(distVicini<=minDist*toll_coeff),"r.");
%     
%     subplot(2,1,3);
%     if(toll_coeff==round(toll_coeff))
%         circle(x(idxTest),y(idxTest),toll_coeff*minDist);
%     end
%     if(sum(distVicini<=minDist*toll_coeff)<=5 && sum(distVicini<=minDist*toll_coeff)>=3)
%         tolprova5=toll_coeff;
%     end
% end
% 
% tolprova=3;
% 
% pastCandSpaz=find(distVicini<=minDist*tolprova);
% subplot(2,1,1);
% hold on;
% for i=pastCandSpaz
%     if(i<=l-samplesPerWin+1)
%         plot(t(i),err(i),"bx");
%     else
%         plot(t(i),min(err(2:end))*2.5,"gs");
%     end
% end
% 
% giroPastCand=-1+zeros(1,length(pastCandSpaz));
% for i=1:length(pastCandSpaz)
%     for j=1:nGiri
%         if(tApproxGiri(j) < t(pastCandSpaz(i)) && t(pastCandSpaz(i))<=tApproxGiri(j+1))
%             giroPastCand(i)=j-1;
%         end
%     end
% end

%% ANALISI 1: QUASI-DISTANZA
%si trova i primi vicini per ogni giro, se ne vede la percentuale giro e
%prende la MEDIA delle percGiro per valutare la percGiro attuale.

%TODO FIX: INVECE DEL PRIMO VICINO, PRENDERE UNA SOVRAPPOSIZIONE DI
%FINESTRE;
%TODO FIX: SERVE DI DISINGUERE I DUE RAMI DELLA CROCE; ESISTONO CASI in cui
%becca come primo vicino il ramo della croce sbagliata!

%NOTA: 
%T=mean(diff(tApproxGiri));
%percGiro=(t(idxTest)-tApproxGiri(end))/T;
%non funziona

percGiro=-1+zeros(1,nGiri);
idxPastCorr=-1*ones(1,nGiri);

figure(fig1);
hold on;
for iGiro=1:nGiri
    %idxPrev=findNearestSpatialCorrespondant(x(idxTest),y(idxTest),iGiro);
    %non funziona: andiamo di closest window
    idxPastCorr(iGiro)=idxPrev;    
    plot(x(idxPrev),y(idxPrev),"rx");
    percGiro(iGiro)=(t(idxPrev)-tApproxGiri(iGiro))/(tApproxGiri(iGiro+1)-tApproxGiri(iGiro));
    if(~(percGiro(iGiro)>=0 &&percGiro(iGiro)<=1))
        disp("inizio")
        findNearestTemporalTimestap(tApproxGiri(iGiro))
        disp("fine")
        findNearestTemporalTimestap(tApproxGiri(iGiro+1))
    end
end
   
percLast=mean(percGiro);
pastCandQuasiDist=zeros(1,nGiri); %se un giro è appena appena completato

for i=1:nGiri%-1
    tCorrispIdeale=tApproxGiri(i)*(1-percLast)+percLast*tApproxGiri(i+1);
    auxvar=abs(t-tCorrispIdeale);
    pastCandQuasiDist(i)= min(find(auxvar==min(auxvar)));
    %va gestito il caso in cuici siano 2 elementi che siano ESATTAMENTE il min
    %può capitare sia quasi 0, occhio agli errori numerici!

end


figure(fig0);
subplot(2,1,1);
hold on;
for i=pastCandQuasiDist
    if(i<=l-samplesPerWin-1)
        plot(t(i),err(i),"kh");
    end
end
figure(fig1);
hold on;
plot(x(pastCandQuasiDist),y(pastCandQuasiDist),"rh");
%% costruzione insieme di punti: 1 quasi-distanza 
pastCand=pastCandQuasiDist;
res1x=x(end-samplesPerWin+1:end);
res1y=y(end-samplesPerWin+1:end);
resQuasiDist=zeros(1,2);
flag=true;

resfig1=figure();
hold on;
title("analisi 1 quasi dist");
plot(res1x,res1y);

for i=pastCand
    lb=max(0,i-samplesPerWin+1);
    rb=min(i+samplesPerWin-1,findNearestTemporalTimestap(tApproxGiri(nGiri)));
    %i left e right bounds sono tra l'inizio di tutto e la fine dell'ultimo
    %giro
    auxvar=x(lb:rb);
    res1x(end+1:end+length(auxvar))=auxvar;
    auxvar=y(lb:rb);
    res1y(end+1:end+length(auxvar))=auxvar;
    if(flag)
        resQuasiDist=auxvar;
        flag=false;
    else
        resQuasiDist(end+1:end+length(auxvar))=auxvar;
    end
    plot(x(lb:rb),y(lb:rb),"x:");
end
plot(x,y,"k.");


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

function [idxPrev]=findClosestPastWindow(idxTest,giro)
%trova il punto più vicono al punto in input
    global x;
    global t;
    global y;
    global tApproxGiri;
    global samplesPerWin;
    inizio=findNearestTemporalTimestap(tApproxGiri(giro));
    fine=findNearestTemporalTimestap(tApproxGiri(giro+1));
    assert(fine-inizio>=samplesPerWin+5);
    
    distlocal=-1*ones(fine-inizio+1-samplesPerWin);
    w0x=x(idxTest-samplesPerWin+1:idxTest);
    w0y=y(idxTest-samplesPerWin+1:idxTest);
    for i=inizio-1+samplesPerWin:length(distlocal)
        wx=x(i-samplesPerWin+1:i);
        wy=y(i-samplesPerWin+1:i);
        %dio cane
        DIO CANE
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