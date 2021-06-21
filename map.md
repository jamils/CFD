# Main

(Will need to store data in JuliaDB)
1. Load constants from input file
2. Get coordinates and build mesh from file
3. Create mesh that includes ghost cells
4. Create temperature distribution
5. Setup variables
## Initialize
6. Set geometry
7. Set primitive variables
8. Set boundary conditions/sources
9. Use MUSCL to get edge states
10. Compute 2d flux
11. Compute residuals
12. Find maxspeed
13. Compute local and global timestep
14. Compute norms and errors for the initial step
15. Get initial residual for normalization
16. Compute error (for MMS case)
17. Output initial arrays
## Runge-Kutta
18. Main time loop - cylcle through time for each cell:
    1.  RK to update interior
    2.  Convert to prim to apply BC
    3.  Set BC
    4.  Use MUSCL scheem with prim
    5.  Compute 2D flux
    6.  Compute residuals
    7.  Compute max speed to get new dt
    8.  Compute timestep local and global with new variables
    9.  Update RK
    10. Compute L2 error
    11. Output timestep arrays
19. Post-process