clear all
close all

H2 = [[0 0 1; -1 0 0; 0 -1 0] [0 -0.06 1.65]'];
H3 = [[0 0 1; -1 0 0; 0 -1 0] [0 0.48 1.65]'];

C1 = [H2; zeros(1,3) 1];
C2 = [H3; zeros(1,3) 1];

figure(1)
R1 = C1(1:3,1:3);
T1 = C1(1:3,end);
R2 = C2(1:3,1:3);
T2 = C2(1:3,end);
quiver3(T1(1),T1(2),T1(3),R1(1,1),R1(2,1),R1(3,1),'r')
hold on
quiver3(T1(1),T1(2),T1(3),R1(1,2),R1(2,2),R1(3,2),'g')
quiver3(T1(1),T1(2),T1(3),R1(1,3),R1(2,3),R1(3,3),'b')


quiver3(T2(1),T2(2),T2(3),R2(1,1),R2(2,1),R2(3,1),'r')
quiver3(T2(1),T2(2),T2(3),R2(1,2),R2(2,2),R2(3,2),'g')
quiver3(T2(1),T2(2),T2(3),R2(1,3),R2(2,3),R2(3,3),'b')

axis equal
grid on
xlabel('x')
ylabel('y')
zlabel('z')