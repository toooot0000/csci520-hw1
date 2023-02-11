<Please submit this file with your solution.>

CSCI 520, Assignment 1

Runan Ye

================

What I have accomplished

1. (Major) implemented the basic jello cube spring structure
2. (Major) implemented dragging interaction. You can drag the cube now!
3. (Minor) implemented Vector3 and matrix-related functions like det()
4. (Minor) implemented ray-casting which is used 

Extra credits:

1. inclined plane collision detection
2. implemented animation interative.
  More about this, check the bottom jello.h where all the variables that are related to dragging are defined. The default drag rate is 1500, the lower the smaller the dragging force is
  Also, the dragging only works with RK4 integrator. My implementation can be very unstable using Euler integrator even without dragging so I didn't bother to test the dragging with it. 
