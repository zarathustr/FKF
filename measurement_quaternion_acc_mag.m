function [q,Jacob]=measurement_quaternion_acc_mag(acc,mag,mag_r,q_)

ax=acc(1);  ay=acc(2);  az=acc(3);
mx=mag(1);  my=mag(2);  mz=mag(3);
mN=mag_r(1);            mD=mag_r(3);

q0=q_(1);   q1=q_(2);   q2=q_(3);   q3=q_(4);

q=zeros(4,1);
Jacob=zeros(4,6);

q(1)= (ay*mD*my + (1 + az)*(1 + mN*mx + mD*mz) + ax*(mD*mx - mN*mz))*q0 + ((mD + az*mD - ax*mN)*my + ay*(1 + mN*mx - mD*mz))*q1 + ...
   (ay*mN*my + ax*(-1 + mN*mx + mD*mz) + (1 + az)*(-(mD*mx) + mN*mz))*q2 + (-((ax*mD + mN + az*mN)*my) + ay*(mD*mx + mN*mz))*q3;

q(2)= ((mD - az*mD - ax*mN)*my + ay*(1 + mN*mx + mD*mz))*q0 + (ay*mD*my - (-1 + az)*(1 + mN*mx - mD*mz) + ax*(mD*mx + mN*mz))*q1 + ...
   ((ax*mD + mN - az*mN)*my + ay*(-(mD*mx) + mN*mz))*q2 + (-(ay*mN*my) + ax*(1 - mN*mx + mD*mz) - (-1 + az)*(mD*mx + mN*mz))*q3;

q(3)= (-(ay*mN*my) - ax*(1 + mN*mx + mD*mz) + (-1 + az)*(mD*mx - mN*mz))*q0 + ((-(ax*mD) + mN - az*mN)*my + ay*(mD*mx + mN*mz))*q1 + ...
   (ay*mD*my + (-1 + az)*(-1 + mN*mx + mD*mz) + ax*(mD*mx - mN*mz))*q2 + ((mD - az*mD + ax*mN)*my + ay*(1 - mN*mx + mD*mz))*q3;

q(4)= ax*(q1 + mN*mx*q1 + mN*my*q2 + mN*mz*q3 + mD*(my*q0 - mz*q1 + mx*q3)) + (1 + az)*(mD*mx*q1 + mD*my*q2 + q3 + mD*mz*q3 - mN*(my*q0 - mz*q1 + mx*q3)) + ... 
   ay*(mN*mz*q0 + mN*my*q1 + q2 - mN*mx*q2 - mD*(mx*q0 + mz*q2 - my*q3));



Jacob(1,1)= -q2 - mN*(mz*q0 + my*q1 - mx*q2) + mD*(mx*q0 + mz*q2 - my*q3);
Jacob(1,2)= q1 + mN*mx*q1 + mN*my*q2 + mN*mz*q3 + mD*(my*q0 - mz*q1 + mx*q3);
Jacob(1,3)= q0 + mN*mx*q0 + mD*mz*q0 + mD*my*q1 - mD*mx*q2 + mN*mz*q2 - mN*my*q3;
Jacob(1,4)= (ax*mD + mN + az*mN)*q0 + ay*mN*q1 + (-((1 + az)*mD) + ax*mN)*q2 + ay*mD*q3;
Jacob(1,5)= ay*mD*q0 + (mD + az*mD - ax*mN)*q1 + ay*mN*q2 - (ax*mD + mN + az*mN)*q3;
Jacob(1,6)= mD*(q0 + az*q0 - ay*q1 + ax*q2) + mN*(-(ax*q0) + q2 + az*q2 + ay*q3);

Jacob(2,1)= q3 - mN*(my*q0 - mz*q1 + mx*q3) + mD*(mx*q1 + my*q2 + mz*q3);
Jacob(2,2)= q0 + mN*mx*q0 + mD*mz*q0 + mD*my*q1 - mD*mx*q2 + mN*mz*q2 - mN*my*q3;
Jacob(2,3)= -((1 + mN*mx)*q1) - mD*(my*q0 - mz*q1 + mx*q3) - mN*(my*q2 + mz*q3);
Jacob(2,4)= ay*(mN*q0 - mD*q2) - (-1 + az)*(mN*q1 + mD*q3) + ax*(mD*q1 - mN*q3);
Jacob(2,5)= mD*(q0 - az*q0 + ay*q1 + ax*q2) - mN*(ax*q0 + (-1 + az)*q2 + ay*q3);
Jacob(2,6)= ay*(mD*q0 + mN*q2) + mD*((-1 + az)*q1 + ax*q3) + mN*(ax*q1 + q3 - az*q3);

Jacob(3,1)= -((1 + mN*mx + mD*mz)*q0) - mD*my*q1 + mD*mx*q2 - mN*mz*q2 + mN*my*q3;
Jacob(3,2)= q3 - mN*(my*q0 - mz*q1 + mx*q3) + mD*(mx*q1 + my*q2 + mz*q3);
Jacob(3,3)= -q2 - mN*(mz*q0 + my*q1 - mx*q2) + mD*(mx*q0 + mz*q2 - my*q3);
Jacob(3,4)= mD*((-1 + az)*q0 + ay*q1 + ax*q2) - mN*(ax*q0 + q2 - az*q2 + ay*q3);
Jacob(3,5)= ay*(-(mN*q0) + mD*q2) - (-1 + az)*(mN*q1 + mD*q3) + ax*(-(mD*q1) + mN*q3);
Jacob(3,6)= mN*(q0 - az*q0 + ay*q1) - ax*(mD*q0 + mN*q2) + mD*((-1 + az)*q2 + ay*q3);

Jacob(4,1)= q1 + mN*mx*q1 + mN*my*q2 + mN*mz*q3 + mD*(my*q0 - mz*q1 + mx*q3);
Jacob(4,2)= q2 + mN*(mz*q0 + my*q1 - mx*q2) - mD*(mx*q0 + mz*q2 - my*q3);
Jacob(4,3)= q3 - mN*(my*q0 - mz*q1 + mx*q3) + mD*(mx*q1 + my*q2 + mz*q3);
Jacob(4,4)= -(ay*(mD*q0 + mN*q2)) + ax*(mN*q1 + mD*q3) + (1 + az)*(mD*q1 - mN*q3);
Jacob(4,5)= (1 + az)*(-(mN*q0) + mD*q2) + ax*(mD*q0 + mN*q2) + ay*(mN*q1 + mD*q3);
Jacob(4,6)= ay*(mN*q0 - mD*q2) + (1 + az)*(mN*q1 + mD*q3) + ax*(-(mD*q1) + mN*q3);

end