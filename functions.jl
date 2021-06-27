function readGrid(imax, jmax)
    file = "griddata.h5";
    grouping = "$imax x $jmax";

    x = h5read(file, grouping*"/x");
    y = h5read(file, grouping*"/y");

    return x, y
end

function cellCenteredGrid(x, y)
    # Returns grid values in the center of the cell

    xc = zeros(jmax, imax);
    yc = zeros(jmax, imax);

    for j ∈ 1:(jmax-1)
        for i ∈ 1:(imax-1)
            xc[j,i] = 0.25 * (x[j,i] + x[j+1,i] + x[j,i+1] + x[j+1,i+1]);
            yc[j,i] = 0.25 * (y[j,i] + y[j+1,i] + y[j,i+1] + y[j+1,i+1]);
        end
    end

    xc[end,end] = 0;
    yc[end,end] = 0;

    return xc, yc
end

function computeArea(xn, yn)
    Ai = zeros(jmax, imax);
    for i ∈ 1:imax
        for j ∈ 1:(jmax-1)
            dx = xn[j,i] - xn[j+1,i];
            dy = yn[j,i] - yn[j+1,i];
            Ai[j,i] = √(dx^2 * dy^2);
        end
    end

    Aj = zeros(jmax, imax);
    for i ∈ 1:(imax-1)
        for j ∈ 1:jmax
            dx = xn[j,i] - xn[j,i+1];
            dy = yn[j,i] - yn[j,i+1];
            Aj[j,i] = √(dx^2 * dy^2);
        end
    end

    return Ai, Aj
end

function computeNormalVectors(xn, yn, Ai, Aj)
    ni_xhat = zeros(jmax, imax);
    ni_yhat = zeros(jmax, imax);
    for i ∈ 1:imax
        for j ∈ 1:(jmax-1)
            ni_xhat[j,i] =  (yn[j+1,i] - yn[j,i]) / Ai[j,i];
            ni_yhat[j,i] = -(xn[j+1,i] - xn[j,i]) / Ai[j,i];
        end
    end

    nj_xhat = zeros(jmax, imax);
    nj_yhat = zeros(jmax, imax);
    for i ∈ 1:(imax-1)
        for j ∈ 1:jmax
            nj_xhat[j,i] = -(yn[j,i+1] - yn[j,i]) / Aj[j,i];
            nj_yhat[j,i] =  (xn[j,i+1] - xn[j,i]) / Aj[j,i];
        end
    end

    return ni_xhat, ni_yhat, nj_xhat, nj_yhat
end

function computeVolume(xn, yn)
    Volume = zeros(jmax, imax);
    for i ∈ 1:(imax-1)
        for j ∈ 1:(jmax-1)
            dx1 = xn[j,i] - xn[j+1,i+1]; # Bottom left to top right
            dy1 = yn[j,i] - yn[j+1,i+1];
            dx2 = xn[j+1,i] - xn[j,i+1]; # Top left to bottom right
            dy2 = yn[j+1,i] = yn[j,i+1];

            Volume[j,i] = 0.5 * abs(dx1*dy2 - dx2*dy1);
        end
    end

    return Volume
end

function extrapCopyCoords(xc, yc)
    xc_g = zeros(jmax+2*num_ghost, imax+2*num_ghost);
    yc_g = zeros(jmax+2*num_ghost, imax+2*num_ghost);
    for i ∈ 1:imax
        for j ∈ 1:jmax
            ig = i + num_ghost;
            jg = j + num_ghost;
            xc_g[jg,ig] = xc[j,i];
            yc_g[jg,ig] = yc[j,i];
        end
    end

    bot = num_ghost + 1; # num_ghost - 1;
    top = num_ghost + jmax;
    left = num_ghost + 1; # num_ghost - 1;
    right = num_ghost + imax;

    # Extrapolate to ghost cells
    for i ∈ num_ghost:(imax + num_ghost)
        for j ∈ 1:num_ghost
            xc_g[bot-j,i] = 2 * xc_g[bot-j+1,i] - xc_g[bot-j+2,i];
            yc_g[bot-j,i] = 2 * yc_g[bot-j+1,i] - yc_g[bot-j+2,i];

            xc_g[top+j,i] = 2 * xc_g[top+j-1,i] - xc_g[top+j-2,i];
            yc_g[top+j,i] = 2 * yc_g[top+j-1,i] - yc_g[top+j-2,i];
        end
    end

    for i ∈ 1:num_ghost
        for j ∈ num_ghost:(jmax + num_ghost)
            xc_g[j,left-i] = 2 * xc_g[j,left-i+1] - xc_g[j,left-i+2]; 
            yc_g[j,left-i] = 2 * yc_g[j,left-i+1] - yc_g[j,left-i+2]; 

            xc_g[j,right+i] = 2 * xc_g[j,right+i-1] - xc_g[j,right+i-2];
            yc_g[j,right+i] = 2 * xc_g[j,right+i-1] - xc_g[j,right+i-2];
        end
    end

    # Set 1 in corners
    for i ∈ 1:num_ghost
        for j ∈ 1:num_ghost
            xc_g[j,i] = 1;
            yc_g[j,i] = 1;
            xc_g[j,right+i] = 1;
            yc_g[j,right+i] = 1;
            xc_g[top+j,right+i] = 1;
            yc_g[top+j,right+i] = 1;
            xc_g[top+j,i] = 1;
            yc_g[top+j,i] = 1;
        end
    end

    return xc_g, yc_g
end

function computeTemperature(V)
    T = zeros(jmax, imax);

    # for i ∈ 1:imax
    #     for j ∈ 1:jmax
    #         T[j,i] = V[j,i,4] / (R*V[j,i,1]);
    #     end
    # end

    @. T = V[:,:,4] / (R*V[:,:,1]);

    return T
end

function prim2cons(V)
    U = zeros(jmaxg, imaxg, 4);

    ρ = V[:,:,1];
    u = V[:,:,2];
    v = V[:,:,3];
    p = V[:,:,4];

    @. ρ = ρ;
    ρu = @. ρ * u;
    ρv = @. ρ * v;
    ρeₜ = @. p / (γ-1) + 0.5 * ρ * u^2 + 0.5 * v^2;

    U[:,:,1] = ρ;
    U[:,:,2] = ρu;
    U[:,:,3] = ρv;
    U[:,:,4] = ρeₜ;

    return U
end

function cons2prim(U)
    V = zeros(jmaxg, imaxg, 4);

    ρ = U[:,:,1];
    ρu = U[:,:,2];
    ρv = U[:,:,3];
    ρeₜ = U[:,:,4];

    @. ρ = ρ;
    @. u = ρu / ρ;
    @. v = ρv / ρ;
    @. p = (γ-1) * (ρeₜ - 0.5 * ρ * u^2 - 0.5 * ρ * v^2); # !!!!!!!!

    V[:,:,1] = ρ;
    V[:,:,2] = u;
    V[:,:,3] = v;
    V[:,:,4] = p;

    return V
end

function computeRes(F, G, S, Ai, Aj, Volume)
    R = zeros(jmax, imax, 4);
    
    for j ∈ 1:jmax
        for i ∈ 1:imax
            @. R[j,i,:] = F[j,i+1,:]*Ai[j,i+1] + G[j+1,i,:]*Aj[j+1,i] - 
                F[j,i,:]*Ai[j,i] - G[j,i,:]*Aj[j,i] - S[j,i,:]*Volume[j,i];
        end
    end

    return R
end

function computeMaxSpeed(V)
    MaxSpeed = zeros(jmax, imax);

    for j ∈ 1:jmax
        for i ∈ 1:imax
            ic = i+num_ghost;
            jc = j+num_ghost;

            MaxSpeed[j,i] = √(V[jc,ic,4]*γ/V[jc,ic,1]);
        end
    end

    return MaxSpeed
end

function computeTimeStep(Volume, Ai, Aj, ni_x, ni_y, nj_x, nj_y, MaxSpeed, dt_min, V)
    dt = zeros(jmax, imax)

    # dt_min = 1e10;

    for j ∈ 1:jmax
        for i ∈ 1:imax
            ic = i+num_ghost;
            jc = j+num_ghost;

            ni_x_avg = 0.5 * (ni_x[j,i] + ni_x[j,i+1]);
            ni_y_avg = 0.5 * (ni_y[j,i] + ni_y[j,i+1]);
            nj_x_avg = 0.5 * (nj_x[j,i] + nj_x[j+1,i]);
            nj_y_avg = 0.5 * (nj_y[j,i] + nj_y[j+1,i]);

            Ai_avg = 0.5 * (Ai[j,i] + Ai[j,i+1]);
            Aj_avg = 0.5 * (Aj[j,i] + Aj[j+1,i]);

            λi = abs(V[jc,ic,2]*ni_x_avg + V[jc,ic,3]*ni_y_avg) + MaxSpeed[j,i];
            λj = abs(V[jc,ic,2]*nj_x_avg + V[jc,ic,3]*nj_y_avg) + MaxSpeed[j,i];

            # Compute local timestep
            dt[j,i] = CFL * Volume[j,i] / (λi * Ai_avg + λj * Aj_avg);

            if dt[j,i] < dt_min
                dt_min = dt[j,i];
            end
        end
    end

    return dt, dt_min
end

function computeL2(Res) # !!!!!!!!
    temp = @. Res^2;
    L2 = √(sum(temp)/4);

    return L2
end

function computeError(V, V_MMS)
    Error = zeros(jmax, imax, 4);

    for j ∈ 1:jmax
        for i ∈ 1:imax
            ic = i+num_ghost;
            jc = j+num_ghost;

            @. Error[j,i,:] = V[jc,ic,:] - V_MMS[jc,ic,:];
        end
    end

    return Error
end

function RungeKutta(U, Res, Volume, dt, dt_min)
    α = zeros(4);
    U_RK = zeros(jmaxg, imaxg, 4);

    if RK_order == 1
        α[1] = 1;
        α[2] = 0;
        α[3] = 0;
        α[4] = 0;

    elseif RK_order == 2
        α[1] = 0.5;
        α[2] = 1;
        α[3] = 0;
        α[4] = 0;

    elseif RK_order == 4
        α[1] = 0.25;
        α[2] = 1/3;
        α[3] = 0.5;
        α[4] = 1;

    end

    for j ∈ jmax
        for i ∈ imax
            if localdt == 0
                dt_val = dt_min;
            else
                dt_val = dt[j,i];
            end

            ic = i+num_ghost;
            jc = j+num_ghost;

            @. U_RK[jc,ic,:] = U[jc,ic,:] - α * dt_val * Res[j,i,:] / Volume[j,i];
        end
    end

    return U_RK
end