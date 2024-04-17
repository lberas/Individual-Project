clear
clc

% Car data (with driver)
% Mass {kg}
m=298;
% Static weight distribution - front
swd_f=0.47;
% Centre of gravity height {m}
cog_h=0.3002;
% Track (average) {m}
t=1.1425;
% Wheelbase {m}
w=1.575;
% Downforce (per corner) {N}
aero=60;

% From tyre data
% Coefficient of friction
cof=1.8;
% Maximum self-aligning torque coefficient {Nm/N}
sat_cof=90/1080;

% Suspension points [x y z] {m}

% Inboard points
% Upper Fore
UF_in=[0.1668 0.2645 0.2563];
% Upper Aft
UA_in=[-0.1899 0.2645 0.2563];
% Lower Fore
LF_in=[0.1333 0.2299 0.1267];
% Lower Aft
LA_in=[-0.1742 0.2299 0.1267];
% Push Rod
PR_in=[0 0.3043 0.595];
% Tie Rod
TR_in=[0.095 0.2061 0.1687];

%Outboard points
% Upper Fore
UF_out=[0 0.4902 0.2998];
% Upper Aft
UA_out=[0 0.4902 0.2998];
% Lower Fore
LF_out=[0 0.4934 0.1372];
% Lower Aft
LA_out=[0 0.4934 0.1372];
% Push Rod
PR_out=[0 0.4536 0.3158];
% Tie Rod
TR_out=[0.065 0.5243 0.1679];

% Member vectors
UF=UF_out-UF_in;
UA=UA_out-UA_in;
LF=LF_out-LF_in;
LA=LA_out-LA_in;
PR=PR_out-PR_in;
TR=TR_out-TR_in;

%Contact patch [x y z]
CP=[0 0.5895 0.0135];

% Maximum braking and cornering acceleration {g}
a_x_max=-cof*(m*9.81+4*aero)/(m*9.81);
a_y_max=-cof*swd_f*(m*9.81+4*aero)/(m*9.81);

% Storing largest forces and largest forces in opposite direction
member_max=ones(6,1);
member_opposite=-ones(6,1);
UBJ_max=ones(3,1);
UBJ_opposite=-ones(3,1);
LBJ_max=ones(3,1);
LBJ_opposite=ones(3,1);
TBJ_max=ones(3,1);
TBJ_opposite=-ones(3,1);


% Entire range of steering angle, negative for inside tyre and positive for
% outside
for s=-40:1:40

    % Rotated unit vectors
    UF_unit=[UF(1)*cosd(s)-UF(2)*sind(s) UF(1)*sind(s)+UF(2)*cosd(s) UF(3)];
    UF_unit=UF_unit/norm(UF_unit);

    UA_unit=[UA(1)*cosd(s)-UA(2)*sind(s) UA(1)*sind(s)+UA(2)*cosd(s) UA(3)];
    UA_unit=UA_unit/norm(UA_unit);

    LF_unit=[LF(1)*cosd(s)-LF(2)*sind(s) LF(1)*sind(s)+LF(2)*cosd(s) LF(3)];
    LF_unit=LF_unit/norm(LF_unit);

    LA_unit=[LA(1)*cosd(s)-LA(2)*sind(s) LA(1)*sind(s)+LA(2)*cosd(s) LA(3)];
    LA_unit=LA_unit/norm(LA_unit);

    PR_unit=[PR(1)*cosd(s)-PR(2)*sind(s) PR(1)*sind(s)+PR(2)*cosd(s) PR(3)];
    PR_unit=PR_unit/norm(PR_unit);

    TR_unit=[TR(1)*cosd(s)-TR(2)*sind(s) TR(1)*sind(s)+TR(2)*cosd(s) TR(3)];
    TR_unit=TR_unit/norm(TR_unit);

    % Range of turning (0 deg) to braking (90 deg); turning only
    % after steering angle >= 10 deg, linear transition
    for phi=0:1:max(90-9*s,0)

        % Forces at the contact patch
        % Downwards force
        Fz=(m*swd_f/2 ...
            -a_x_max*swd_f*cog_h*m*(sind(phi)*cosd(s)+cosd(phi)*sind(s))/...
            (2*w)-a_x_max*(1-swd_f)*cog_h*m*sind(phi)/(2*w)...
            -a_y_max*cog_h*m*(-sind(phi)*sind(s)+cosd(phi)*cosd(s))/...
            (2*t))*9.81+aero;
        % Braking will be limited by inside tyre loading
        Fz_inner=(m*swd_f/2 ...
            -a_x_max*swd_f*cog_h*m*(sind(phi)*cosd(-abs(s))...
            -cosd(phi)*sind(-abs(s)))/(2*w)...
            -a_x_max*(1-swd_f)*cog_h*m*sind(phi)/(2*w)...
            -a_y_max*cog_h*m*(sind(phi)*sind(-abs(s))...
            +cosd(phi)*cosd(-abs(s)))/(2*t))*9.81+aero;
        % Braking force
        Fx=-Fz_inner*cof*sind(phi);
        % Turning force
        Fy=-Fz*cof*cosd(phi);

        % Moments to be resolved by matrix equation, SAT depends on
        % steering angle magintude and direction
        if s==0 % No SAT at zero steering angle
            moments=cross(CP,[Fx Fy Fz])+[0 0 0];
        elseif abs(s)<4 % SAT reaches maximum at 4 deg steering angle
            moments=cross(CP,[Fx Fy Fz])+[0 0 sat_cof*Fz*s/4];
        else % SAT maintains maximum, s/abs(s) for direction
            moments=cross(CP,[Fx Fy Fz])+[0 0 sat_cof*Fz*s/abs(s)];
        end

        % Matrix equation setup
        A=[UF_unit' UA_unit' LF_unit' LA_unit' PR_unit' TR_unit';...
            cross(UF_out,UF_unit)' cross(UA_out,UA_unit)'...
            cross(LF_out,LF_unit)' cross(LA_out,LA_unit)'...
            cross(PR_out,PR_unit)' cross(TR_out,TR_unit)'];

        % Matrix equation solution and storing largest values
        members=A\[Fx Fy Fz moments]';
        for i = 1:length(members)
            if abs(members(i))>abs(member_max(i))
                if members(i)/member_max(i)<0
                    member_opposite(i)=member_max(i);
                end
                member_max(i)=members(i);
            end
        end

        % Finding forces at ball joints and storing largest values
        UBJ=sum(-[UF_unit; UA_unit; PR_unit].*members([1 2 5]),1);
        for i=1:length(UBJ)
            if abs(UBJ(i))>abs(UBJ_max(i))
                if UBJ(i)/UBJ_max(i)<0
                    UBJ_opposite(i)=UBJ_max(i);
                end
                UBJ_max(i)=UBJ(i);
            end
        end

        LBJ=sum(-[LF_unit; LA_unit].*members([3 4]),1);
        for i=1:length(LBJ)
            if abs(LBJ(i))>abs(LBJ_max(i))
                if LBJ(i)/LBJ_max(i)<0
                    LBJ_opposite(i)=LBJ_max(i);
                end
                LBJ_max(i)=LBJ(i);
            end
        end

        TBJ=sum(-TR_unit.*members(6),1);
        for i=1:length(TBJ)
            if abs(TBJ(i))>abs(TBJ_max(i))
                if TBJ(i)/TBJ_max(i)<0
                    TBJ_opposite(i)=TBJ_max(i);
                end
                TBJ_max(i)=TBJ(i);
            end
        end
    end
end

member_max
UBJ_max
LBJ_max
TBJ_max
member_opposite
UBJ_opposite
LBJ_opposite
TBJ_opposite

