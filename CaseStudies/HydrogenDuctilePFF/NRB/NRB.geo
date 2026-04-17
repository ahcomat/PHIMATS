// Gmsh project created on Sun Feb 08 10:31:55 2026
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {5e-3, 0, 0, 1.0};
//+
Point(3) = {5e-3, 6.282e-3, 0, 1.0};
//+
Point(4) = {3e-3, 15e-3, 0, 1.0};
//+
Point(5) = {3e-3, 35e-3, 0, 1.0};
//+
Point(6) = {5e-3, 43.718e-3, 0, 1.0};
//+
Point(7) = {5e-3, 50e-3, 0, 1.0};
//+
Point(8) = {0, 50e-3, 0, 1.0};
//+
Point(9) = {0, 25e-3, 0, 1.0};
//+
Point(10) = {23e-3, 15e-3, 0, 1.0};
//+
Point(11) = {23e-3, 35e-3, 0, 1.0};
//+
Point(13) = {3e-3, 27.236e-3, 0, 1.0};
//+
Point(14) = {3e-3, 22.764e-3, 0, 1.0};
//+
Point(15) = {5e-3, 25e-3, 0, 1.0};
//+
Circle(1) = {14, 15, 13};
//+
Circle(2) = {3, 10, 4};
//+
Circle(3) = {5, 11, 6};
//+
Line(4) = {1, 2};
//+
Line(5) = {2, 3};
//+
Line(6) = {4, 14};
//+
Line(7) = {13, 5};
//+
Line(8) = {6, 7};
//+
Line(9) = {7, 8};
//+
Line(10) = {8, 1};
//+
Line(11) = {9, 15};
//+
Curve Loop(1) = {7, 3, 8, 9, 10, 4, 5, 2, 6, 1};
//+
Plane Surface(1) = {1};
//+
Physical Surface("NRB", 12) = {1};
//+
Physical Curve("Symmetry_Axis", 13) = {10};
//+
Physical Curve("Bottom", 14) = {4};
//+
Physical Curve("Top", 15) = {9};
//+
Physical Point("Fix", 16) = {1};

// --- Distance field fillet ---
Field[1] = Distance;
Field[1].EdgesList = {1,11}; 

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 1e-5;     // fine mesh near notch
Field[2].SizeMax = 1e-3;     // coarser far away
Field[2].DistMin = 0.0;      // start refining right at the notch
Field[2].DistMax = 2e-2;     // blending distance 

// --- Combine all refinement fields ---
Field[10] = Min;
Field[10].FieldsList = {2};  // combine all thresholds

// --- Use this field as background mesh size ---
Background Field = 10;

Recombine Surface {1};
