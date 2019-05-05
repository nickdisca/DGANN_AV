L = 4;// channel length
H = 1; // height
x0 = 1/6.0; // starting point of ramp
h=1.0/200.0;

Point(1) = {0,0,0,h};
Point(2) = {x0,0,0,h};
Point(3) = {L,0,0,h};
Point(4) = {L,H,0,h};
Point(5) = {0,H,0,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,1};

Line Loop(1) = {1,2,3,4,5};
Plane Surface(1) = {1};


Physical Surface(100) = {1};

Physical Line(101) = {5}; // left
Physical Line(102) = {3}; // right
Physical Line(103) = {1}; // bottom before ramp
Physical Line(104) = {2}; // ramp
Physical Line(105) = {4}; // top


