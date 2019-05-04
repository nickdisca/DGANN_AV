L = 1;// channel length
H = 1; // height

n = 100;  // number of points on each boundary

Point(1) = {0,0,0};
Point(2) = {L,0,0};
Point(3) = {L,H,0};
Point(4) = {0,H,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,3,4};

Transfinite Line{1,-3} = n;
Transfinite Line{2,-4} = n;

Physical Surface(100000) = {1};

Physical Line(100001) = {4}; // left
Physical Line(100002) = {2}; // right
Physical Line(100003) = {1}; // bottom
Physical Line(100004) = {3}; // top

//Periodic lines
Periodic Line {1} = {-3};
Periodic Line {2} = {-4};
