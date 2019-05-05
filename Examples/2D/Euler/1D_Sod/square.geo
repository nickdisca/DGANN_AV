L = 6;// channel length
H = 0.2; // height
h=0.04;

Point(1) = {-L/2,-H/2,0,h};
Point(2) = {L/2,-H/2,0,h};
Point(3) = {L/2,H/2,0,h};
Point(4) = {-L/2,H/2,0,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};


Physical Surface(100) = {1};

Physical Line(101) = {4}; // left
Physical Line(102) = {2}; // right
Physical Line(103) = {1}; // bottom
Physical Line(104) = {3}; // top

//Periodic lines
Periodic Line {1} = {-3};

