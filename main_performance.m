clear all;
close all;
clc;

addpath('quaternion_library');
addpath('linspecer');

data=load('flat_motion3.txt');
conv=diag([-1 1 -1]);
Accelerometer=data(:,4:6)*conv;
Magnetometer=data(:,10:12)*conv;
Gyroscope=data(:,7:9)*180/pi*conv;
euler_true=data(:,1:3)*conv;

Sigma_g=0.001*eye(3);
Sigma_a=0.01*eye(3);
Sigma_m=0.01*eye(3);
Pk=0.001*eye(4);
q=[1;0;0;0];

dt=1/500;
len=length(data(:,1));
time=dt*(1:len);

quaternion=zeros(len,4);
quaternion_measurement=zeros(len,4);
STD=zeros(len,4);

for i=1:2000
    q0=q(1);
    q1=q(2);
    q2=q(3);
    q3=q(4);
    
    Accelerometer(1,:)=Accelerometer(1,:)./norm(Accelerometer(1,:));
    Magnetometer(1,:)=Magnetometer(1,:)./norm(Magnetometer(1,:));
    
    mD=dot(Accelerometer(1,:),Magnetometer(1,:));
    mN=sqrt(1-mD^2);
    
    [q, Jacob]=measurement_quaternion_acc_mag(Accelerometer(1,:),Magnetometer(1,:),[mN,0,mD], q);
    q=q./norm(q);
end

for i=1:len
    q0=q(1);
    q1=q(2);
    q2=q(3);
    q3=q(4);
    
    wx=Gyroscope(i,1)*pi/180;
    wy=Gyroscope(i,2)*pi/180;
    wz=Gyroscope(i,3)*pi/180;
    
    Accelerometer(i,:)=Accelerometer(i,:)./norm(Accelerometer(i,:));
    Magnetometer(i,:)=Magnetometer(i,:)./norm(Magnetometer(i,:));
    
    mD=dot(Accelerometer(i,:),Magnetometer(i,:));
    mN=sqrt(1-mD^2);
    
    omega4=[0,-wx,-wy,-wz;
        wx,0,wz,-wy;
        wy,-wz,0,wx;
        wz,wy,-wx,0];
   
    Phi=eye(4)+dt/2*omega4;
    
    Dk=[q1 q2 q3;
        -q0 -q3 -q2;
        q2 -q0 -q1;
        -q2 q1 -q0];
    Xi=dt*dt/4*Dk*Sigma_g*Dk';
    
    [qy, Jacob]=measurement_quaternion_acc_mag(Accelerometer(i,:),Magnetometer(i,:),[mN,0,mD], q);
    qy=qy./norm(qy);
    
    Eps=Jacob*[Sigma_a,zeros(3,3);zeros(3,3),Sigma_m]*Jacob';
    
    std=sqrt([Eps(1,1),Eps(2,2),Eps(3,3),Eps(4,4)]);
    
    q_=q;
    Pk_ = Pk;
    [q , Pk] = kalman_update(q_, qy, Pk_, Phi, Xi, Eps);
    
    
    q=q./norm(q);
    quaternion(i,:)=q';
    quaternion_measurement(i,:)=qy';
    STD(i,:)=std;
    
end


euler=quatern2euler(quaternConj(quaternion))*180/pi;

for i=1:len
    if euler(i,3)<0
        euler(i,3)=euler(i,3)+360;
    end
end

figure(1);
subplot(3,1,1);
plot(time,euler_true(:,1),'LineWidth',2); hold on;
plot(time,euler(:,1),'LineWidth',1); hold off;
xlabel('Time (s)');
ylabel('Roll (deg)');
legend('Reference','Estimated');

subplot(3,1,2);
plot(time,euler_true(:,2),'LineWidth',2); hold on;
plot(time,euler(:,2),'LineWidth',1); hold off;
xlabel('Time (s)');
ylabel('Pitch (deg)');
legend('Reference','Estimated');

subplot(3,1,3);
plot(time,euler_true(:,3)+180,'LineWidth',2); hold on;
plot(time,euler(:,3),'LineWidth',1); hold off;
xlabel('Time (s)');
ylabel('Yaw (deg)');
legend('Reference','Estimated');

Colors=linspecer(2);
figure(2);
subplot(4,1,1);
plot(time,quaternion_measurement(:,1),'LineWidth',2,'Color',Colors(2,:)); hold on;
plot(time,quaternion_measurement(:,1)+3*STD(:,1),'LineWidth',1,'Color',Colors(1,:)); hold on;
plot(time,quaternion_measurement(:,1)-3*STD(:,1),'LineWidth',1,'Color',Colors(1,:)); hold off;
xlabel('Time (s)');
ylabel('{q_0}');
legend('Measurement Quaternion','3\sigma Bound');

subplot(4,1,2);
plot(time,quaternion_measurement(:,2),'LineWidth',2,'Color',Colors(2,:)); hold on;
plot(time,quaternion_measurement(:,2)+3*STD(:,2),'LineWidth',1,'Color',Colors(1,:)); hold on;
plot(time,quaternion_measurement(:,2)-3*STD(:,2),'LineWidth',1,'Color',Colors(1,:)); hold off;
xlabel('Time (s)');
ylabel('{q_1}');
legend('Measurement Quaternion','3\sigma Bound');

subplot(4,1,3);
plot(time,quaternion_measurement(:,3),'LineWidth',2,'Color',Colors(2,:)); hold on;
plot(time,quaternion_measurement(:,3)+3*STD(:,3),'LineWidth',1,'Color',Colors(1,:)); hold on;
plot(time,quaternion_measurement(:,3)-3*STD(:,3),'LineWidth',1,'Color',Colors(1,:)); hold off;
xlabel('Time (s)');
ylabel('{q_2}');
legend('Measurement Quaternion','3\sigma Bound');

subplot(4,1,4);
plot(time,quaternion_measurement(:,4),'LineWidth',2,'Color',Colors(2,:)); hold on;
plot(time,quaternion_measurement(:,4)+3*STD(:,4),'LineWidth',1,'Color',Colors(1,:)); hold on;
plot(time,quaternion_measurement(:,4)-3*STD(:,4),'LineWidth',1,'Color',Colors(1,:)); hold off;
xlabel('Time (s)');
ylabel('{q_3}');
legend('Measurement Quaternion','3\sigma Bound');


figure(3);
subplot(3,1,1);
plot(time,data(:,7:9)*conv,'LineWidth',1);
legend('x','y','z');
xlabel('Time (s)');
ylabel('Angular Rate (rad/s)');

subplot(3,1,2);
plot(time,data(:,4:6)*conv,'LineWidth',1);
legend('x','y','z');
xlabel('Time (s)');
ylabel('Acceleration (G)');

subplot(3,1,3);
plot(time,data(:,10:12)*conv,'LineWidth',1);
legend('x','y','z');
xlabel('Time (s)');
ylabel('Magnetic Field (Gauss)');