close all
plot3(0,0,0);
hold on;
n = size(pos ,1 );
%ref = ones(n)*260;
%t=times(1:n);
%plot(t,ref, '--','MarkerSize',0.5,'MarkerFaceColor','b');
for i=[1:round(n/100):n]
    
%     plot3(pos(i,1),pos(i,2),pos(i,3),'or','MarkerSize',2,'MarkerFaceColor','r')
%     plot3(tgtPos(i,1),tgtPos(i,2),tgtPos(i,3),'or','MarkerSize',2,'MarkerFaceColor','b')
     plot3(pos(i,1),pos(i,2),pos(i,3),"rx");
     plot3(tgtPos(i,1),tgtPos(i,2),tgtPos(i,3),"+k");

    %plot3(ArrayRM1(i),ArrayRM2(i),'or','MarkerSize',2,'MarkerFaceColor','b')
    %axis([0 400 70 170])
%     axis([0 400 100 135])
%     set(gca,'FontSize',12)
%     xlabel('Distance (s)')
%     ylabel('Height (deg)')
%     %legend('Target','Missile')
    pause(0.05);
    i
end
