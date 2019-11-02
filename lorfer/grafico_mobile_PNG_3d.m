close all
plot3(0,0,0);
hold on;
n = size(ArrayT ,2 );
%ref = ones(n)*260;
%t=times(1:n);
%plot(t,ref, '--','MarkerSize',0.5,'MarkerFaceColor','b');
for i=[1:58]
    
%     plot3(pos_dronei,1),pos_dronei,2),pos_dronei,3),'or','MarkerSize',2,'MarkerFaceColor','r')
%     plot3(pos_target(i,1),pos_target(i,2),pos_target(i,3),'or','MarkerSize',2,'MarkerFaceColor','b')
     plot3(array_pos_drone(i,1),array_pos_drone(i,2),array_pos_drone(i,3),"rx");
     plot3(array_pos_target(i,1),array_pos_target(i,2),array_pos_target(i,3),"+k");

    %plot3(ArrayRM1(i),ArrayRM2(i),'or','MarkerSize',2,'MarkerFaceColor','b')
    %axis([0 400 70 170])
%     axis([0 400 100 135])
%     set(gca,'FontSize',12)
%     xlabel('Distance (s)')
%     ylabel('Height (deg)')
%     %legend('Target','Missile')
    pause(0.03);
end