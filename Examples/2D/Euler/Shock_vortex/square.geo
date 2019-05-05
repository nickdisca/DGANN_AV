h=0.01;

Point(1) = {0,0,0,h};
Point(2) = {2,0,0,h};
Point(3) = {2,1,0,h};
Point(4) = {0,1,0,h};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};


Physical Surface(100) = {1};

Physical Line(101) = {4}; // inflow
Physical Line(102) = {1}; // bottom
Physical Line(103) = {2}; // outflow
Physical Line(104) = {3}; // top


