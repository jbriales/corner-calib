function result=PluckerCoordinates2PluckerMatrix(L)

result=[0 L(6) -L(5) L(1); -L(6) 0 L(4) L(2); L(5) -L(4) 0 L(3); -L(1) -L(2) -L(3) 0];
