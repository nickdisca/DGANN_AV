L = 10;// channel length
H = 10; // height

n = 80;  // number of points on each boundary

Point(1) = {-L/2,-H/2,0};
Point(2) = {L/2,-H/2,0};
Point(3) = {L/2,H/2,0};
Point(4) = {-L/2,H/2,0};

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