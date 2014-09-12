close all;
clc;
clear all;

% P = [1 -1 -1.5; 1 0.5 2; 3 4 2];
P = rand(3,3);
P(1:2,:) = P(1:2,:) - 0.5;
P(3,:) = P(3,:) + 2; 
R = (pi/6)*(rand(3,1)-0.5);
R = angle2dcm(R(3),R(2),R(1));
t = rand(3,1)-0.5;
T = [R t; 0 0 0 1];
Ti = inv(T);

x = Ti*[P; ones(1,3)];
x = x(1:3,:)./(ones(3,1)*[norm(x(1:3,1))...
                          norm(x(1:3,2))...
                          norm(x(1:3,3))]);

[Te, s1, s2, s3] = p3p(x,P,'grunert');

figure;
scatter3(P(1,:),P(2,:),P(3,:),'k');
hold on;
%PlotAxes(eye(4),gca,1,0);
plot3([Ti(1,4) P(1,1)],[Ti(2,4) P(2,1)],[Ti(3,4) P(3,1)],'k');
plot3([Ti(1,4) P(1,2)],[Ti(2,4) P(2,2)],[Ti(3,4) P(3,2)],'k');
plot3([Ti(1,4) P(1,3)],[Ti(2,4) P(2,3)],[Ti(3,4) P(3,3)],'k');


for i=1:size(Te,3)
    Tei(:,:,i) = inv(Te(:,:,i));
    Pe = Tei(:,:,i)*([x.*(ones(3,1)*[s1(i) s2(i) s3(i)]);ones(1,3)]);
 
    hold on;
    PlotAxes(Tei(:,:,i),gca,0.5,0);
    plot3([Te(1,4,i) Pe(1,1)],[Te(2,4,i) Pe(2,1)],[Te(3,4,i) Pe(3,1)],'r');
    plot3([Te(1,4,i) Pe(1,2)],[Te(2,4,i) Pe(2,2)],[Te(3,4,i) Pe(3,2)],'r');
    plot3([Te(1,4,i) Pe(1,3)],[Te(2,4,i) Pe(2,3)],[Te(3,4,i) Pe(3,3)],'r');
end

axis equal



