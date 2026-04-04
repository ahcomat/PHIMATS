SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {63.5e-3, 0, 0, 1.0};
//+
Point(3) = {63.5e-3, 0, 0, 1.0};
//+
Point(4) = {63.5e-3, 20.5e-3, 0, 1.0};
//+
Point(5) = {50.8e-3, 20.5e-3, 0, 1.0};
//+
Point(6) = {50.8e-3, 28.5e-3, 0, 1.0};
//+
Point(7) = {50.8e-3, 32.5e-3, 0, 1.0};
//+
Point(8) = {50.8e-3, 40.5e-3, 0, 1.0};
//+
Point(9) = {63.8e-3, 40.5e-3, 0, 1.0};
//+
Point(10) = {63.8e-3, 61e-3, 0, 1.0};
//+
Point(11) = {0, 61e-3, 0, 1.0};
//+
Point(12) = {34.5e-3, 28.5e-3, 0, 1.0};
//+
Point(13) = {34.5e-3, 32.5e-3, 0, 1.0};
//+
Point(14) = {30.4e-3, 30.53e-3, 0, 1.0};
//+
Point(15) = {50.8e-3, 11.5e-3, 0, 1.0};
//+
Point(16) = {44.45e-3, 11.5e-3, 0, 1.0};
//+
Point(17) = {57.15e-3, 11.5e-3, 0, 1.0};
//+
Point(18) = {50.8e-3, 49.5e-3, 0, 1.0};
//+
Point(19) = {44.45e-3, 49.5e-3, 0, 1.0};
//+
Point(20) = {57.15e-3, 49.5e-3, 0, 1.0};
//+
Point(21) = {22.859e-3, 30.5e-3, 0, 1e-7};
//+
Point(22) = {21.0e-3, 30.5e-3, 0, 1.0};
//+
Point(23) = {31.5e-3, 30.5e-3, 0, 1.0};
//+
Point(24) = {30.4e-3, 30.58e-3, 0, 1.0};
//+
Point(25) = {25.025e-3, 30.502e-3, 0, 1.0};
//+
Point(26) = {25.025e-3, 30.4986e-3, 0, 1.0};

//+
Circle(1) = {17, 15, 16};
//+
Circle(2) = {16, 15, 17};
//+
Circle(3) = {19, 18, 20};
//+
Circle(4) = {20, 18, 19};
//+
Line(5) = {1, 2};
//+
Line(6) = {2, 4};
//+
Line(7) = {4, 5};
//+
Line(8) = {5, 6};
//+
Line(9) = {6, 12};
//+
Line(10) = {12, 14};

Spline(11) = {25, 21, 26};
//+
Line(13) = {24, 13};
//+
Line(14) = {13, 7};
//+
Line(15) = {7, 8};
//+
Line(16) = {8, 9};
//+
Line(17) = {9, 10};
//+
Line(18) = {10, 11};
//+
Line(19) = {11, 1};
//+
Line(21) = {21, 22};
//+
Line(22) = {24, 25};
//+
Line(23) = {26, 14};

//+
Curve Loop(1) = {5, 6, 7, 8, 9, 10, 23, 11, 22, 13, 14, 15, 16, 17, 18, 19};
//+
Curve Loop(2) = {2, 1};
//+
Curve Loop(3) = {4, 3};
//+
Plane Surface(1) = {1, 2, 3};

//+
Physical Surface("CT", 21) = {1};
//+
Physical Point("Fix", 22) = {2};
//+
Physical Point("Upper_gauge", 23) = {7};
//+
Physical Point("Lower_gauge", 24) = {6};
//+
Physical Curve("Lower_curve", 25) = {1};
//+
Physical Curve("Upper_curve", 26) = {3};
//+
Physical Curve("Surface", 27) = {5, 6, 1, 2, 7, 8, 9, 10, 23, 11, 22, 13, 14, 15, 16, 4, 3, 17, 18, 19};

// --- Distance field center ---
Field[1] = Distance;
Field[1].EdgesList = {21}; 

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 1e-5;    // fine mesh near notch
Field[2].SizeMax = 2e-3;    // coarser far away
Field[2].DistMin = 0;       // start refining right at the notch
Field[2].DistMax = 2.5e-2;  // blending distance 

// --- Distance field fillet ---
Field[3] = Distance;
Field[3].EdgesList = {11}; 

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = 1e-5;     // fine mesh near edge
Field[4].SizeMax = 2e-3;     // coarser far away
Field[4].DistMin = 0.0;      // start refining right at the edge
Field[4].DistMax = 2e-2;     // blending distance 

// --- Distance field fillet ---
Field[5] = Distance;
Field[5].EdgesList = {1,2,3,4}; 

Field[6] = Threshold;
Field[6].InField = 5;
Field[6].SizeMin = 0.001;    // fine mesh near edge
Field[6].SizeMax = 2e-3;     // coarser far away
Field[6].DistMin = 0.0;      // start refining right at the edge
Field[6].DistMax = 1e-2;     // blending distance 

// --- Combine all refinement fields ---
Field[10] = Min;
Field[10].FieldsList = {2,4,6};  // combine all thresholds

// --- Use this field as background mesh size ---
Background Field = 10;

Recombine Surface {1};

