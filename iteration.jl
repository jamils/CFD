for k ∈ 1:RK_order
    if n == 1 && k == 1
        dt_min = 1e10;
        # Set primitive variables

        println("Initialize")

        V = initialize(xc_g, yc_g);

        T = computeTemperature(V);

        if case == 1 # MMS Curvilinear mesh
            S = solveSourceMMS(xc, yc);
            V_MMS = solveSolutionMMS(xc_g, yc_g);
            V = setBCMMS(V_MMS);

            # !!!!!!!! Output to DB

        else
            S = zeros(jmax, imax);
            V = setBC(ni_x, ni_y, nj_x, nj_y, T);
        end

        U = prim2cons(V);
        
        VL, VR, VB, VT = MUSCL(V);

        F, G = compute2DFlux(ni_x, ni_y, nj_x, nj_y, VL, VR, VB, VT);

        Res = computeRes(F, G, S, Ai, Aj, Volume); # Residuals

        MaxSpeed = computeMaxSpeed(V);

        dt, dt_min = computeTimeStep(
            Volume, Ai, Aj, ni_x, ni_y, nj_x,nj_y, MaxSpeed, dt_min, V);
        
        U_RK = RungeKutta(U, Res, Volume, dt, dt_min);
        V = cons2prim(U_RK);
        U = U_RK;

        L2 = computeL2(Res);

        if case == 1 # Compute error for MMS case
            Error = computeError(V, V_MMS);
            Inter = computeL2(Error);
        end
            println("Starting Main Loop")
    end

    println("RK = $k")

    println("RK")
    U_RK = RungeKutta(U, Res, Volume, dt, dt_min);
    V = cons2prim(U_RK);
    println("Temp")
    T = computeTemperature(V);
    if case != 1
        println("BC")
        V = setBC(ni_x, ni_y, nj_x, nj_y, T);
    end

    println("MUSCL")
    VL, VR, VB, VT = MUSCL(V);
    println("2DFlux")
    F, G = compute2DFlux(ni_x, ni_y, nj_x, nj_y, VL, VR, VB, VT);
    println("Residuals")
    Res = computeRes(F, G, S, Ai, Aj, Volume);
end

MaxSpeed = computeMaxSpeed(V);
println("dt")
dt, dt_min = computeTimeStep(
    Volume, Ai, Aj, ni_x, ni_y, nj_x, nj_y, MaxSpeed, dt_min, V);

U = U_RK;
println("L2")
L2 = computeL2(Res);

if case == 1
    Error = computeError(V, V_MMS);
    Iter = computeL2(Error);
end