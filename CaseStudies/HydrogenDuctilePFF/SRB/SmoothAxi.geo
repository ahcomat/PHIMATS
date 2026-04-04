// Gmsh project created on Sun Feb 15 05:13:10 2026
SetFactory("OpenCASCADE");
// Dimensions
R_gauge = 0.003;   
R_grip  = 0.006;    
L_total = 0.05;      // Total length of the bar
f_rad   = 0.003;     // Fillet radius (Difference between radii)

// Vertical positions
Y_mid   = L_total / 2.0;    // Point of necking (0.02)
Y_gauge_start = Y_mid - 0.015; 
Y_gauge_end   = Y_mid + 0.015;

// Points - Starting from Bottom Left (0,0)
Point(1) = {0, 0, 0, 0.001};                   // Bottom on Y-axis
Point(2) = {R_grip, 0, 0, 0.001};              // Bottom Right (Grip)
Point(3) = {R_grip, Y_gauge_start - f_rad, 0, 0.001}; // Start of lower fillet

// Lower Fillet
Point(4) = {R_grip, Y_gauge_start, 0, 0.0005}; // Center for lower fillet
Point(5) = {R_gauge, Y_gauge_start, 0, 0.0005}; // End of lower fillet

// Necking Zone
Point(6) = {R_gauge, Y_mid, 0, 0.0001}; // Necking point at mid-height

// Upper Fillet
Point(7) = {R_gauge, Y_gauge_end, 0, 0.0005};   // Start of upper fillet
Point(8) = {R_grip, Y_gauge_end, 0, 0.0005};    // Center for upper fillet
Point(9) = {R_grip, Y_gauge_end + f_rad, 0, 0.001}; // End of upper fillet

Point(10) = {R_grip, L_total, 0, 0.001};       // Top Right
Point(11) = {0, L_total, 0, 0.001};            // Top Left (on Y-axis)

Point(12) = {0, Y_mid, 0, 0.0001}; // Opposite to Necking point at mid-height

// Lines and Curves
Line(1) = {1, 2};                              // Bottom edge
Line(2) = {2, 3};                              // Lower grip edge
Circle(3) = {3, 4, 5};                         // Lower Fillet (3=start, 4=center, 5=end)
Line(4) = {5, 6};                              // Gauge section with necking
Line(10) = {6, 7};                             // Gauge section with necking
Circle(5) = {7, 8, 9};                         // Upper Fillet (7=start, 8=center, 9=end)
Line(6) = {9, 10};                             // Upper grip edge
Line(7) = {10, 11};                            // Top edge
Line(8) = {11, 1};                             // Axis of Symmetry (Y-axis)
//+
Line(9) = {12, 6};

// Surface definition
Curve Loop(1) = {1, 2, 3, 4, 10, 5, 6, 7, 8};
Plane Surface(1) = {1};

// Boundary Conditions
Physical Curve("Bottom") = {1};
Physical Curve("Top") = {7};
Physical Curve("Symmetry_Axis") = {8};
Physical Surface("SmoothAxi") = {1};
Physical Point("Fix", 10) = {1};

Physical Point("Lower_gauge", 11) = {5};
Physical Point("Upper_gauge", 12) = {7};
Physical Curve("Surface", 13) = {5, 4, 10, 3, 6, 2};

// --- Distance field center ---
Field[1] = Distance;
Field[1].EdgesList = {9}; 

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.5e-5;   // fine mesh near notch
Field[2].SizeMax = 1e-3;     // coarser far away
Field[2].DistMin = 0.0;      // start refining right at the notch
Field[2].DistMax = 2e-2;     // blending distance 

// --- Distance field fillet ---
Field[3] = Distance;
Field[3].EdgesList = {5,3}; 

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = 1e-4;     // fine mesh near notch
Field[4].SizeMax = 1e-3;     // coarser far away
Field[4].DistMin = 0.0;      // start refining right at the notch
Field[4].DistMax = 2e-2;     // blending distance 

// --- Distance field fillet ---
Field[5] = Distance;
Field[5].EdgesList = {4, 10}; 

Field[6] = Threshold;
Field[6].InField = 5;
Field[6].SizeMin = 1e-4;     // fine mesh near notch
Field[6].SizeMax = 1e-3;     // coarser far away
Field[6].DistMin = 0.0;      // start refining right at the notch
Field[6].DistMax = 2e-2;     // blending distance 

// --- Combine all refinement fields ---
Field[10] = Min;
Field[10].FieldsList = {2,4,6};  // combine all thresholds

// --- Use this field as background mesh size ---
Background Field = 10;

Recombine Surface {1};

