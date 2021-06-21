function MUSCL(V)
    VL, VR = computeUpwindVLR(V);
    VB, VT = computeUpwindVBT(V);

    return VL, VR, VB, VT
end

function applyLimiter(r⁺, r⁻)
    ψ⁺ = zeros(4);
    ψ⁻ = zeros(4);

    if limiter == 0 # no limiter
        @. ψ⁺ = 1;
        @. ψ⁻ = 1;

    elseif limiter == 1 # vanLeer
        @. ψ⁺ = (r⁺ + abs(r⁺)) / (1+abs(r⁺));
        @. ψ⁻ = (r⁻ + abs(r⁻)) / (1+abs(r⁻));

    elseif limiter == 2 # vanAlbada
        @. ψ⁺ = (r⁺ + r⁺^2) / (1 + r⁺^2);
        @. ψ⁻ = (r⁻ + r⁻^2) / (1 + r⁻^2);

    elseif limiter == 3 # Ospre
        @. ψ⁺ = 1.5 * (r⁺^2 + r⁺) / (1 + r⁺ + r⁺^2);
        @. ψ⁻ = 1.5 * (r⁻^2 + r⁻) / (1 + r⁻ + r⁻^2);

    elseif limiter == 4 # monotized center least diffusive
        @. ψ⁺ = maximum(0, minimum(2.0 * r⁺, minimum(0.5*(1+r⁺),2)));
        @. ψ⁻ = maximum(0, minimum(2.0 * r⁻, minimum(0.5*(1+r⁻),2)));

    elseif limiter == 5 # minmod
        @. ψ⁺ = maximum(0, minimum(1, r⁺));
        @. ψ⁻ = maximum(0, minimum(1, r⁻));

    elseif limiter == 6 # superbee
        @. ψ⁺ = maximum(0, maximum(minimum(2*r⁺, 1), minimum(r⁺, 2)));
        @. ψ⁻ = maximum(0, maximum(minimum(2*r⁻, 1), minimum(r⁻, 2)));

    end

    return ψ⁺, ψ⁻
end

function computeUpwindVLR(V)
    ni_g = imax;
    nj_g = jmax;
    ni = imax - 2*num_ghost;
    nj = jmax - 2*num_ghost;

    Ψ⁺ = zeros(nj, ni_g+1, 4);
    Ψ⁻ = zeros(nj, ni_g+1, 4);

    δ = 1e-6;

    for j ∈ 1:nj
        for i ∈ (num_ghost-2):(ni+num_ghost)
            j_cells = num_ghost + j;
            den = @. V[j_cells,i+1,:] - V[j_cells,i,:];
            for i ∈ 1:4
                if den[i] < δ || den[i] == 0
                    den[i] = δ;
                end
            end
            # Calc slopes for both positive and negative direction
            r⁺ = @. (V[j_cells,i+2,:] - V[j_cells,i+1,:])/den;
            r⁻ = @. (V[j_cells,i,:] - V[j_cells,i-1,:])/den;

            # Applying limiters
            Ψ⁺[j,i+1,:], Ψ⁻[j,i+1,:] = applyLimiter(r⁺, r⁻);
        end
    end

    VL = zeros(jmax, imax, 4);
    VR = zeros(jmax, imax, 4);

    # Fill in left and right states
    for j ∈ 1:jmax
        for i ∈ 1:imax
            i_int = num_ghost + i;
            i_cells = (num_ghost-1) + i;
            j_cells = num_ghost + j;

            @. VL[j,i,:] = V[j_cells,i_cells,:] +
                0.25 * ε * ((1-κ)*Ψ⁺[j,i_int-1,:] * 
                (V[j_cells,i_cells,:] - V[j_cells,i_cells-1,:])) + 
                (1+κ) * Ψ⁻[j,i_int,:] * 
                (V[j_cells,i_cells+1,:] - V[j_cells,i_cells,:]);

            @. VR[j,i,:] = V[j_cells,i_cells+1,:] - 
                0.25 * ε * ((1-κ)*Ψ⁻[j,i_int+1,:] *
                (V[j_cells,i_cells+2,:] - V[j_cells,i_cells+1,:])) + 
                (1+κ) * Ψ⁺[j,i_int,:] * 
                (V[j_cells,i_cells+1,:] - V[j_cells,i_cells,:]);
        end
    end

    return VL, VR
end

function computeUpwindVBT(V)
    ni_g = imax;
    nj_g = jmax;
    ni = imax - 2*num_ghost;
    nj = jmax - 2*num_ghost;

    Ψ⁺ = zeros(nj, ni_g+1, 4);
    Ψ⁻ = zeros(nj, ni_g+1, 4);

    δ = 1e-6;

    for i ∈ 1:ni
        for j ∈ (num_ghost-2):(nj+num_ghost)
            i_cells = num_ghost + i;
            den = @. V[j+1,i_cells,:] - V[j,i_cells,:];
            for i ∈ 1:4
                if den[i] < δ || den[i] == 0
                    den[i] = δ;
                end
            end
            # Calc slopes for both positive and negative direction
            r⁺ = @. (V[j+2,i_cells,:] - V[j+1,i_cells,:])/den;
            r⁻ = @. (V[j,i_cells,:] - V[j-1,i_cells,:])/den;

            # Applying limiters
            Ψ⁺[j+1,i,:], Ψ⁻[j+1,i,:] = applyLimiter(r⁺, r⁻);
        end
    end

    VB = zeros(jmax, imax, 4);
    VT = zeros(jmax, imax, 4);

    # Fill in left and right states
    for i ∈ 1:imax
        for j ∈ 1:jmax
            j_int = num_ghost + j;
            j_cells = (num_ghost-1) + j;
            i_cells = num_ghost + i;

            @. VB[j,i,:] = V[j_cells,i_cells,:] +
                0.25 * ε * ((1-κ)*Ψ⁺[j_int-1,i,:] * 
                (V[j_cells,i_cells,:] - V[j_cells-1,i_cells,:])) + 
                (1+κ) * Ψ⁻[j_int,i,:] * 
                (V[j_cells+1,i_cells,:] - V[j_cells,i_cells,:]);

            @. VT[j,i,:] = V[j_cells+1,i_cells,:] - 
                0.25 * ε * ((1-κ)*Ψ⁻[j_int+1,i,:] *
                (V[j_cells+2,i_cells,:] - V[j_cells+1,i_cells,:])) + 
                (1+κ) * Ψ⁺[j_int,i,:] * 
                (V[j_cells+1,i_cells,:] - V[j_cells,i_cells,:]);
        end
    end

    return VB, VT
end