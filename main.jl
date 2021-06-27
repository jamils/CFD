include("includes.jl")

const imax = 53;
const jmax = 17;
const num_ghost = 3; # 3
const CFL = 0.6;
const nmax = 1e6;
const RK_order = 4;
const R = 287.0;
const γ = 1.4;

const localdt = 0;
const case = 1;
const supersonic = 0;
const limiter = 0;
const upwind = 1;
const κ = 0.25;
const ε = 0.5;

const imaxg = imax + 2*num_ghost;
const jmaxg = jmax + 2*num_ghost;

x, y = readGrid(imax, jmax); # !!!!! HEAVY MODIFY
xc, yc = cellCenteredGrid(x, y);

xc_g, yc_g = extrapCopyCoords(xc, yc);

# Initialize the domain
# Set geometry

xn = x;
yn = y;

Ai, Aj = computeArea(xn, yn);

ni_x, ni_y, nj_x, nj_y = computeNormalVectors(xn, yn, Ai, Aj);

Volume = computeVolume(xn, yn);

# !!!!!! OUTPUT TO DB

# Set primitive variables

println("Initialize")

V = initialize(xc_g, yc_g);

T = computeTemperature(V);

if case == 1 # MMS Curvilinear mesh
    S = solveSourceMMS(xc, yc);
    V_MMS = solveSolutionMMS(xc_g, yc_g);
    V= setBCMMS(V_MMS);

    # !!!!!!!! Output to DB

else
    S = zeros(jmax, imax);

    V = setBC(ni_x, ni_y, nj_x, nj_y, T);
end

U = prim2cons(V);
U_RK = R;

VL, VR, VB, VT = MUSCL(V);

F, G = compute2DFlux(ni_x, ni_y, nj_x, nj_y, VL, VR, VB, VT);

Res = computeRes(F, G, S, Ai, Aj, Volume); # Residuals

MaxSpeed = computeMaxSpeed(V);

dt, dt_min = computeTimeStep(
    Volume, Ai, Aj, ni_x, ni_y, nj_x,nj_y, MaxSpeed, dt_min, V);

L2 = computeL2(Res);

if case == 1 # Compute error for MMS case
    Error = computeError(V, V_MMS);
    Inter = computeL2(Error);
end

# !!!!!! Use function to output initial data to DB

println("Starting main loop")

# RungeKutta Iteration
for n ∈ 1:nmax
    for k ∈ 1:RK_order
        U_RK = RungeKutta(U, Res, Volume, dt, dt_min);
        V = cons2prim(U_RK);
        T = computeTemperature(V);
        if case != 1
            V = setBC(ni_x, ni_y, nj_x, nj_y, T);
        end

        VL, VR, VB, VT = MUSCL(V);
        F, G = compute2DFlux(ni_x, ni_y, nj_x, nj_y, VL, VR, VB, VT);
        Res = computeRes(F, G, S, Ai, Aj, Volume);
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

    # !!!!!!!! OUTPUT TO DB MORE
end

println("Done!!!")