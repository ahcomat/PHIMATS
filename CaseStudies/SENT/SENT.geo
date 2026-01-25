// Gmsh project created on Wed Sep 24 17:46:30 2025
SetFactory("OpenCASCADE");

t = 8.5e-3;
a = 0.27*t;
L = 52e-3;
r = 0.2e-3; 

//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {t, 0, 0, 1.0};
//+
Point(3) = {t, (L/2)-r, 0, 1.0};
//+
Point(4) = {t-0.27*t, (L/2)-r, 0, 1.0};
//+
Point(5) = {t-0.27*t, (L/2), 0, 1.0};
//+
Point(6) = {t-0.27*t, (L/2)+r, 0, 1.0};
//+
Point(7) = {t, (L/2)+r, 0, 1.0};
//+
Point(8) = {t, L, 0, 1.0};
//+
Point(9) = {0, L, 0, 1.0};
//+
Point(10) = {0, (L/2)+5*r, 0, 1.0};
//+
Point(11) = {0, (L/2)-5*r, 0, 1.0};


//+
Circle(1) = {4, 5, 6};
//+
Line(2) = {11, 10};
//+
Line(3) = {4, 3};
//+
Line(4) = {6, 7};
//+
Line(5) = {11, 1};
//+
Line(6) = {1, 2};
//+
Line(7) = {2, 3};
//+
Line(8) = {7, 8};
//+
Line(9) = {8, 9};
//+
Line(10) = {9, 10};
//+
Curve Loop(1) = {8, 9, 10, -2, 5, 6, 7, -3, 1, 4};
//+
Curve Loop(2) = {8, 9, 10, -2, 5, 6, 7, -3, 1, 4};
//+
Plane Surface(1) = {2};
//+
Physical Surface("SENT", 11) = {1};

Field[1] = Distance;
Field[1].EdgesList = {1,2}; 

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 1e-5;     // fine mesh at curves
Field[2].SizeMax = 1e-3;     // coarser far away
Field[2].DistMin = 0.0;      // start refining right at curves
Field[2].DistMax = 2e-2;     // blending distance 

Field[3] = Distance;
Field[3].EdgesList = {4,3};

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = 5e-5;    
Field[4].SizeMax = 1e-3;
Field[4].DistMin = 0.0;
Field[4].DistMax = 1e-2;

Field[10] = Min;
Field[10].FieldsList = {2, 4};  // combine all thresholds

Background Field = 10;
//+
Physical Point("corner", 12) = {1};
//+
Physical Curve("bottom", 13) = {6};
//+
Physical Curve("top", 14) = {9};
