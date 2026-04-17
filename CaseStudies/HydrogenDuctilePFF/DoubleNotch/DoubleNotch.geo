// Gmsh project created on Sun Mar 22 10:00:10 2026
SetFactory("OpenCASCADE");
//+
Point(1) = {1e-3, 1e-3, 0, 1.0};
//+
Point(2) = {3.291e-3, 0, 0, 1.0};
//+
Point(3) = {0, 3.291e-3, 0, 1.0};
//+
Point(4) = {10e-3, 0, 0, 1.0};
//+
Point(5) = {0, 3.291e-3, 0, 1.0};
//+
Point(6) = {9.0e-3, 9.0e-3, 0, 1.0};
//+
Point(7) = {10e-3, 7.268e-3, 0, 1.0};
//+
Point(8) = {7.268e-3, 10e-3, 0, 1.0};
//+
Point(10) = {0, 10e-3, 0, 1.0};
//+
Point(11) = {0.0027675802767761, 0.0027675802767761, 0, 1.0};
//+
Point(12) = {0.00758581, 0.00758581, 0, 1.0};

//+
Circle(1) = {2, 1, 11};
//+
Circle(2) = {12, 6, 8};
//+
Circle(3) = {11, 1, 3};
//+
Circle(4) = {7, 6, 12};
//+
Line(5) = {2, 4};
//+
Line(6) = {4, 7};
//+
Line(7) = {8, 10};
//+
Line(8) = {10, 3};
//+
Curve Loop(1) = {6, 4, 2, 7, 8, -3, -1, 5};
//+
Plane Surface(1) = {1};
//+
Physical Surface("DoubleNotch", 9) = {1};
//+
Physical Curve("Left", 10) = {8};
//+
Physical Curve("Top", 11) = {7};
//+
Physical Curve("Bottom", 12) = {5};
//+
Physical Curve("Right", 13) = {6};
//+

//+
Line(9) = {8, 2};


// --- Distance field fillet ---
Field[1] = Distance;
Field[1].EdgesList = {1,2,9}; 

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 3e-5;     // fine mesh near notch
Field[2].SizeMax = 2e-3;     // coarser far away
Field[2].DistMin = 0.0;      // start refining right at the notch
Field[2].DistMax = 2e-2;     // blending distance 

// --- Combine all refinement fields ---
Field[10] = Min;
Field[10].FieldsList = {2};  // combine all thresholds

// --- Use this field as background mesh size ---
Background Field = 10;

Recombine Surface {1};
