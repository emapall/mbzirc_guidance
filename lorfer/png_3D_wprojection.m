clear all;
close all;
n=0; BETA=-0.; v_drone_start = 30.; v_target_start = 10.; XNT = .966; 
HEDEG = -30.; XNP = 5.; T=0; S=0;
pos_start_drone = [ 0, 100., 0 ]; pos_start_target = [ 400, 100, 10];

v_target = [ v_target_start*cos(BETA), v_target_start*sin(BETA), 0 ];
v_drone = [ 5, 0, 0 ];
pos_drone = pos_start_drone;
pos_target = pos_start_target;

los = pos_target - pos_drone; %line of sight 
d = norm(los);

%Projection on XY plane and YZ
plane_xy = [ 0, 0, 1];
plane_yz = [1, 0, 0];
los_xy = los - (los*plane_xy')*plane_xy;
d_xy = norm(los_xy);
lambda_xy = atan2(los(1), los(2));
los_yz = los - (los*plane_yz')*plane_yz;
d_yz = norm(los_yz);
lambda_yz = atan2(los(2), los(3));
v_TD = v_target - v_drone;
vc_xy = -(los_xy*v_TD')/d_xy; %closing velocity on XY
vc_yz = -(los_yz*v_TD')/d_yz;

while d>=1
    if d < 100
        H=.0002;
    else
        H=.01;
    end
    BETAOLD=BETA;
    pos_targetOLD = pos_target;
    pos_droneOLD = pos_drone;
    v_droneOLD = v_drone;
    STEP = 1;
    FLAG = 0;
    while STEP <=1
        if FLAG==1
                STEP=2;
                BETA=BETA+H*BETAD;
                pos_target = pos_target + v_target*H;
                pos_drone = pos_drone + v_drone*H;
                v_drone = v_drone + a_drone*H;
                T=T+H;
        end
            los = pos_target - pos_drone; %line of sight
            d = norm(los);
            
            los_xy = los - (los*plane_xy')*plane_xy;
            d_xy = norm(los_xy);
            lambda_xy = atan2(los(1), los(2));
            
            los_yz = los - (los*plane_yz')*plane_yz;
            d_yz = norm(los_yz);
            lambda_yz = atan2(los(2), los(3));
            
            v_TD = v_target - v_drone;
            vc_xy = -(los_xy*v_TD')/d_xy;
            vc_yz = -(los_yz*v_TD')/d_yz;
            
            lambda_xy_dot=(los(1)*v_TD(2)-los(2)*v_TD(1))/(d_xy^2);
            lambda_yz_dot=(los(2)*v_TD(3)-los(3)*v_TD(2))/(d_yz^2);
            acc_perp_xy = XNP*vc_xy*lambda_xy_dot;
            acc_perp_yz = XNP*vc_yz*lambda_yz_dot;
            %a_xy = [-acc_perp_xy*sin(lambda_xy), acc_perp_xy*cos(lambda_xy), 0];
            a_yz = [0, -acc_perp_yz*sin(lambda_yz), acc_perp_yz*cos(lambda_yz)];
            a_xy=0;
            a_drone = a_xy + a_yz;
            v_target = [ v_target_start*cos(BETA), v_target_start*sin(BETA), 0 ];
            BETAD=XNT/v_target_start;
        FLAG=1;
    end
    FLAG=0;
    BETA=.5*(BETAOLD+BETA+H*BETAD);
    pos_target = 0.5*(pos_targetOLD+pos_target+H*v_target);
    pos_drone = 0.5*(pos_droneOLD+pos_drone+H*v_drone);
    v_drone = 0.5 * (v_droneOLD+v_drone+H*a_drone);
    S=S+H;
    if S >=.09999
            S=0.;
        n=n+1;
        ArrayT(n)=T;
        array_pos_drone(n,:) = pos_drone;
        array_pos_target(n,:) = pos_target;
        array_v_drone(n,:) = v_drone;
        array_vc_yz(n) = vc_yz;
        array_lambda_yz_dot(n) = lambda_yz_dot;
        array_acc_yz(n,:) = a_yz;
        array_acc_xy(n,:) = a_xy;
        %ArrayXNCG(n)=XNC/32.2;
        array_d(n)=d;
        array_acc_drone(n,:) = a_drone;
    end
end

