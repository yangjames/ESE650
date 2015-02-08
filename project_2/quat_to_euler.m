function [phi,theta,psi] = quat_to_euler(q)
phi = atan2(2*(q(1)*q(2)+q(3)*q(4)),1-2*(q(3)^2+q(4)^2));
theta = asin(2*(q(1)*q(3)-q(4)*q(2)));
psi = atan2(2*(q(1)*q(4)+q(2)*q(3)),1-2*(q(3)^2+q(4)^2));