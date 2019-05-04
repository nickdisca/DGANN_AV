L = 3;// channel length
H = 3; // height
h=0.02;

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


Physical Surface(100000) = {1};

Physical Line(100001) = {4}; // left
Physical Line(100002) = {2}; // right
Physical Line(100003) = {1}; // bottom
Physical Line(100004) = {3}; // top

//Periodic lines
Periodic Line {1} = {-3};
Periodic Line {2} = {-4};

