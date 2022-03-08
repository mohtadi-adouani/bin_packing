\ENCODING=ISO-8859-1
\Problem name: TP1

Minimize
 obj1: y_i0 + y_i1 + y_i2 + y_i3 + y_i4 + y_i5 + y_i6
Subject To
 bin_sac_i0:    - 10 y_i0 + 3 x_i0_j0 + 4 x_i0_j1 + 4 x_i0_j2 + 3 x_i0_j3
                + 3 x_i0_j4 + 2 x_i0_j5 + x_i0_j6 <= 0
 bin_sac_i1:    - 10 y_i1 + 3 x_i1_j0 + 4 x_i1_j1 + 4 x_i1_j2 + 3 x_i1_j3
                + 3 x_i1_j4 + 2 x_i1_j5 + x_i1_j6 <= 0
 bin_sac_i2:    - 10 y_i2 + 3 x_i2_j0 + 4 x_i2_j1 + 4 x_i2_j2 + 3 x_i2_j3
                + 3 x_i2_j4 + 2 x_i2_j5 + x_i2_j6 <= 0
 bin_sac_i3:    - 10 y_i3 + 3 x_i3_j0 + 4 x_i3_j1 + 4 x_i3_j2 + 3 x_i3_j3
                + 3 x_i3_j4 + 2 x_i3_j5 + x_i3_j6 <= 0
 bin_sac_i4:    - 10 y_i4 + 3 x_i4_j0 + 4 x_i4_j1 + 4 x_i4_j2 + 3 x_i4_j3
                + 3 x_i4_j4 + 2 x_i4_j5 + x_i4_j6 <= 0
 bin_sac_i5:    - 10 y_i5 + 3 x_i5_j0 + 4 x_i5_j1 + 4 x_i5_j2 + 3 x_i5_j3
                + 3 x_i5_j4 + 2 x_i5_j5 + x_i5_j6 <= 0
 bin_sac_i6:    - 10 y_i6 + 3 x_i6_j0 + 4 x_i6_j1 + 4 x_i6_j2 + 3 x_i6_j3
                + 3 x_i6_j4 + 2 x_i6_j5 + x_i6_j6 <= 0
 objectUsed_i0: x_i0_j0 + x_i1_j0 + x_i2_j0 + x_i3_j0 + x_i4_j0 + x_i5_j0
                + x_i6_j0  = 1
 objectUsed_i1: x_i0_j1 + x_i1_j1 + x_i2_j1 + x_i3_j1 + x_i4_j1 + x_i5_j1
                + x_i6_j1  = 1
 objectUsed_i2: x_i0_j2 + x_i1_j2 + x_i2_j2 + x_i3_j2 + x_i4_j2 + x_i5_j2
                + x_i6_j2  = 1
 objectUsed_i3: x_i0_j3 + x_i1_j3 + x_i2_j3 + x_i3_j3 + x_i4_j3 + x_i5_j3
                + x_i6_j3  = 1
 objectUsed_i4: x_i0_j4 + x_i1_j4 + x_i2_j4 + x_i3_j4 + x_i4_j4 + x_i5_j4
                + x_i6_j4  = 1
 objectUsed_i5: x_i0_j5 + x_i1_j5 + x_i2_j5 + x_i3_j5 + x_i4_j5 + x_i5_j5
                + x_i6_j5  = 1
 objectUsed_i6: x_i0_j6 + x_i1_j6 + x_i2_j6 + x_i3_j6 + x_i4_j6 + x_i5_j6
                + x_i6_j6  = 1
Bounds
 0 <= y_i0 <= 1
 0 <= y_i1 <= 1
 0 <= y_i2 <= 1
 0 <= y_i3 <= 1
 0 <= y_i4 <= 1
 0 <= y_i5 <= 1
 0 <= y_i6 <= 1
 0 <= x_i0_j0 <= 1
 0 <= x_i0_j1 <= 1
 0 <= x_i0_j2 <= 1
 0 <= x_i0_j3 <= 1
 0 <= x_i0_j4 <= 1
 0 <= x_i0_j5 <= 1
 0 <= x_i0_j6 <= 1
 0 <= x_i1_j0 <= 1
 0 <= x_i1_j1 <= 1
 0 <= x_i1_j2 <= 1
 0 <= x_i1_j3 <= 1
 0 <= x_i1_j4 <= 1
 0 <= x_i1_j5 <= 1
 0 <= x_i1_j6 <= 1
 0 <= x_i2_j0 <= 1
 0 <= x_i2_j1 <= 1
 0 <= x_i2_j2 <= 1
 0 <= x_i2_j3 <= 1
 0 <= x_i2_j4 <= 1
 0 <= x_i2_j5 <= 1
 0 <= x_i2_j6 <= 1
 0 <= x_i3_j0 <= 1
 0 <= x_i3_j1 <= 1
 0 <= x_i3_j2 <= 1
 0 <= x_i3_j3 <= 1
 0 <= x_i3_j4 <= 1
 0 <= x_i3_j5 <= 1
 0 <= x_i3_j6 <= 1
 0 <= x_i4_j0 <= 1
 0 <= x_i4_j1 <= 1
 0 <= x_i4_j2 <= 1
 0 <= x_i4_j3 <= 1
 0 <= x_i4_j4 <= 1
 0 <= x_i4_j5 <= 1
 0 <= x_i4_j6 <= 1
 0 <= x_i5_j0 <= 1
 0 <= x_i5_j1 <= 1
 0 <= x_i5_j2 <= 1
 0 <= x_i5_j3 <= 1
 0 <= x_i5_j4 <= 1
 0 <= x_i5_j5 <= 1
 0 <= x_i5_j6 <= 1
 0 <= x_i6_j0 <= 1
 0 <= x_i6_j1 <= 1
 0 <= x_i6_j2 <= 1
 0 <= x_i6_j3 <= 1
 0 <= x_i6_j4 <= 1
 0 <= x_i6_j5 <= 1
 0 <= x_i6_j6 <= 1
Binaries
 y_i0  y_i1  y_i2  y_i3  y_i4  y_i5  y_i6  x_i0_j0  x_i0_j1  x_i0_j2  x_i0_j3 
 x_i0_j4  x_i0_j5  x_i0_j6  x_i1_j0  x_i1_j1  x_i1_j2  x_i1_j3  x_i1_j4 
 x_i1_j5  x_i1_j6  x_i2_j0  x_i2_j1  x_i2_j2  x_i2_j3  x_i2_j4  x_i2_j5 
 x_i2_j6  x_i3_j0  x_i3_j1  x_i3_j2  x_i3_j3  x_i3_j4  x_i3_j5  x_i3_j6 
 x_i4_j0  x_i4_j1  x_i4_j2  x_i4_j3  x_i4_j4  x_i4_j5  x_i4_j6  x_i5_j0 
 x_i5_j1  x_i5_j2  x_i5_j3  x_i5_j4  x_i5_j5  x_i5_j6  x_i6_j0  x_i6_j1 
 x_i6_j2  x_i6_j3  x_i6_j4  x_i6_j5  x_i6_j6 
End
