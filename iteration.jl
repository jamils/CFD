for k ∈ 1:RK_order
    println("n = $n, k = $k")
    global U, U_RK, Res, Volume, dt, dt_min, V, T, VL, VR, VB, VT, F, G, MaxSpeed, L2, Error, Iter;

    U_RK = RungeKutta(U, Res, Volume, dt, dt_min);
    V = cons2prim(U_RK);
    T = computeTemperature(V);
    if case != 1
        V = setBC(ni_x, ni_y, nj_x, nj_y, T);
    end

    VL, VR, VB, VT = MUSCL(V);
    F, G = compute2DFlux(ni_x, ni_y, nj_x, nj_y, VL, VR, VB, VT);
    Res = computeRes(F, G, S, Ai, Aj, Volume);
    U = prim2cons(V);
end

MaxSpeed = computeMaxSpeed(V);
dt, dt_min = computeTimeStep(
    Volume, Ai, Aj, ni_x, ni_y, nj_x, nj_y, MaxSpeed, dt_min, V);

U = U_RK;
L2 = computeL2(Res);

if case == 1
    Error = computeError(V, V_MMS);
    Iter = computeL2(Error);
end

n += 1;

print(" ")
heatmap(U[:,:,2])
# lineplot(V[:,begin,4])

# !!!!!!!! OUTPUT TO DB MORE