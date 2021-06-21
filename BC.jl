function setBC(ni_x, ni_y, nj_x, nj_y, T)
    ni = imax - 2*num_ghost;
    nj = jmax - 2*num_ghost;

    if case == 1
        error("Wrong BC!")
        
    elseif case == 2
        joint = 10 * 2^(mesh);

        Lower_Begin = [0, 0];
        Lower_End = [0, ni-1];

        Upper_Begin = [nj-1, joint];
        Upper_End = [nj-1, ni-1];

        Inlet_Begin = [nj-1, 0];
        Inlet_End = [nj-1, joint-1];

        Outlet_Begin = [0, ni];
        Outlet_End = [nj-1, ni];

        Sym_Begin = [0,0];
        Sym_End = [nj-1, 0];

        V = symmetricBC(V, Sym_Begin, Sym_End);
        V = inletBC(V, Inlet_Begin, Inlet_End);
        V = outletBC(V, Outlet_Begin, Outlet_End);
        V = slipwallBC(V, Upper_Begin, Upper_End, ni_x, ni_y, nj_x, nj_y, T);
        V = slipwallBC(V, Lower_Begin, Lower_End, ni_x, ni_y, nj_x, nj_y, T);

    elseif case == 3 # NACA Airfoil
        iBegin = 4 * 2^mesh;
        iEnd = 20 * 2^mesh - 1;

        # !!!!! May need to increase values (Or just 0s)

        Farfield1_Begin = [0, 0];
        Farfield1_End = [nj-1, 0];
        
        Farfield2_Begin = [nj-1, 0];
        Farfield2_End = [nj-1, ni-1];

        Farfield3_Begin = [0, ni-1];
        Farfield3_End = [nj-1, ni-1];

        Periodic_A_Begin = [0, 0];
        Periodic_A_End = [0, iBegin-1];

        Periodic_B_Begin = [0, ni-1];
        Periodic_B_End = [0, iEnd];

        Airfoil_Begin = [0, iBegin];
        Airfoil_End = [0, iEnd];

        V = farfieldBC(V, Farfield1_Begin, Farfield1_End);
        V = farfieldBC(V, Farfield2_Begin, Farfield2_End);
        V = farfieldBC(V, Farfield3_Begin, Farfield3_End);

        V = slipwallBC(V, Airfoil_Begin, Airfoil_End, ni_x, ni_y, nj_x, nj_y, T);
        V = periodicBC(V, Periodic_A_Begin, Periodic_A_End, Periodic_B_Begin, Periodic_B_End);

    end

    return V
end

function symmetricBC(V, B, E)
    iB = B[2] + num_ghost;
    jB = B[1] + num_ghost;

    for j ∈ 1:(E[1]-B[1])
        for i ∈ 1:num_ghost
            I = iB - 1;
            J = jB;

            V[J+j,I-i,1] = V[J+j,I+i+1,1];
            V[J+j,I-i,2] = V[J+j,I+i+1,2];
            V[J+j,I-i,3] = V[J+j,I+i+1,3];
            V[J+j,I-i,4] = V[J+j,I+i+1,4];
        end
    end

    return V
end

function inletBC(V, B, E)
    iB = B[2] + num_ghost;
    jB = B[1] + num_ghost;

    M = 4;
    p₀ = 12270; # [Pa]
    T₀ = 217; # [K]
    ρ₀ = p₀ / (R*T₀); # [kg/m³]
    u₀ = M * √(γ*R*T₀);
    v₀ = 0;

    for j ∈ 1:num_ghost
        for i ∈ 1:(E[2]-B[2])
            I = iB;
            J = jB + 1;

            V[J+j,I+i,1] = ρ₀;
            V[J+j,I+i,2] = u₀;
            V[J+j,I+i,3] = v₀;
            V[J+j,I+i,4] = p₀;
        end
    end

    return V
end

function outletBC(V, B, E)
    iB = B[2] + num_ghost;
    jB = B[1] + num_ghost;

    for j ∈ 1:(E[1]-B[1])
        for i ∈ 1:num_ghost
            I = iB;
            J = jB;

            @. V[J+j,I+i, :] = 2 * V[J+j,I+i-1,:] - V[J+j,I+i-2,:];
        end
    end

    return V
end

function slipwallBC(V, B, E, ni_x, ni_y, nj_x, nj_y, T)
    ni_g = imax;
    nj_g = jmax;
    ni = ni_g - 2*num_ghost;
    nj = nj_g - 2*num_ghost;

    iB = B[2] + num_ghost;
    jB = B[1] + num_ghost;
    iE = E[2] + num_ghost;
    jE = E[1] + num_ghost;

    if iB == iE
        if iB == num_ghost
            sign = -1;
            BC_index = B[2];

        elseif iB == (ni + num_ghost-1)
            sign = 1;
            BC_index = B[2] + 1;

        end

        for j ∈ 1:(E[1]-B[1])
            for i ∈ 1:num_ghost
                I = iB + sign;
                J = jB;
                nx = ni_x[B[1]+j, BC_index];
                ny = ni_y[B[1]+j, BC_index];
                u = V[J+j,iB-sign*i,2];
                v = V[J+j,iB-sign*i,3];
                Temp = T[J+j,iB-sign*i];

                V[J+j,I+sign*i,2] = - nx*(u*nx+v*ny) - ny*(-u*ny+v*nx);
                V[J+j,I+sign*i,3] = - ny*(u*nx+v*ny) + nx*(-u*ny+v*nx);

                V[J+j,I+sign*i,4] = 2 * V[J+j,I+sign*(i-1),4] - V[J+j,I+sign*(i-2),4];
                T[J+j,I+sign*i] = Temp;
                V[J+j,I+sign*i,1] = V[J+j,I+sign*i,4] / (R*T[J+j,I+sign*i]);
            end
        end

    elseif jB == jE
        if jB == num_ghost
            sign = -1;
            BC_index = B[1];

        elseif jB == (nj+num_ghost-1)
            sign = 1;
            BC_index = B[1] + 1;

        end

        for j ∈ 1:num_ghost
            for i ∈ 1:(E[2]-B[2])
                I = iB;
                J = jB;
                nx = nj_x[BC_index,B[2]+i];
                ny = nj_y[BC_index,B[2]+i];

                u = V[jB-sign*j,I+i,2];
                v = V[jB-sign*j,I+i,3];
                Temp = T[jB-sign*j,I+i];

                V[J+sign*j,I+i,2] = -nx*(u*nx+v*ny) - ny*(-u*ny+v*nx);
                V[J+sign*j,I+i,3] = -ny*(u*nx+v*ny) + nx*(-u*ny+v*nx);
                V[J+sign*j,I+i,4] = 2*V[J+sign*(j-1),I+i,4] - V[J+sign*(j-2),I+i,4];
                
                T[J+sign*j,I+i] = Temp;
                V[J+sign*j,I+i,1] = V[J+sign*j,I+i,4] / (R*T[J+sign*j,I+i]);
            end
        end
    end

    return V
end

function farfieldBC(V, B, E)
    ni_g = imax;
    nj_g = jmax;
    ni = ni_g - 2*num_ghost;
    nj = nj_g - 2*num_ghost;

    iB = B[2] + num_ghost;
    jB = B[1] + num_ghost;
    iE = E[2] + num_ghost;
    jE = E[1] + num_ghost;

    if AngleOfAttack == 1
        M = 0.84;
        p₀ = 65855.8; # [Pa]
        T₀ = 300; # [K]
        AOA = 0; # [ᵒ]

    elseif AngleOfAttack == 2
        M = 0.75;
        p₀ = 67243.5; # [Pa]
        T₀ = 300; # [K]
        AOA = 8; # [ᵒ]

    end

    ρ₀ = P₀ / (R*T₀);
    u₀ = M*√(γ*R*T₀)*cos(AOA*π/180);
    v₀ = M*√(γ*R*T₀)*sin(AOA*π/180);

    if iB == iE
        if iB == num_ghost
            sign = -1;
        elseif iB == (ni+num_ghost-1)
            sign = 1;
        end

        for j ∈ 1:(E[1]-B[1])
            for i ∈ 1:num_ghost
                I = iB + sign;
                J = jB;

                V[J+j,I+sign*i,1] = ρ₀;
                V[J+j,I+sign*i,2] = u₀;
                V[J+j,I+sign*i,3] = v₀;
                V[J+j,I+sign*i,4] = p₀;
            end
        end

    elseif jB == jE
        if jB == num_ghost
            sign = -1;
        elseif jB == (nj+num_ghost-1)
            sign = 1;
        end

        for j ∈ 1:num_ghost
            for i ∈ 1:(E[2]-B[2])
                I = iB;
                J = jB + sign;

                V[J+sign*j,I+i,1] = ρ₀;
                V[J+sign*j,I+i,2] = u₀;
                V[J+sign*j,I+i,3] = v₀;
                V[J+sign*j,I+i,4] = p₀;
            end
        end
    end

    return V
end

function periodicBC(V, BA, EA, BB, EB)
    ni_g = imax;
    nj_g = jmax;
    ni = ni_g - 2*num_ghost;
    nj = nj_g - 2*num_ghost;

    iBA = BA[2] + num_ghost;
    jBA = BA[1] + num_ghost;
    iEA = EA[2] + num_ghost;
    jEA = EA[1] + num_ghost;

    iBB = BB[2] + num_ghost;
    jBB = BB[1] + num_ghost;
    iEB = EB[2] + num_ghost;
    jEB = EB[1] + num_ghost;

    if iBA == iEA
        if iBA == num_ghost
            sign = -1;
        elseif iBA == (ni+num_ghost-1)
            sign = 1;
        end

    elseif jBA == jEA
        L = abs(iBA - iEA) + 1;

        if jBA == num_ghost
            sign = -1;
        elseif jBA == (nj+num_ghost-1)
            sign = 1;
        end

        if iEB > iBB
            Bdir = 1;
        elseif iEB < iBB
            Bdir = -1;
        end

        for j ∈ 1:num_ghost
            for i ∈ 1:L
                IA = iBA;
                IB = iBB;

                JA = JBA + sign;
                JB = JBB + sign;

                @. V[JA+sign*j,IA+i,:] = V[jBB-sign*j,IB+Bdir*i,:];
                @. V[JB+sign*j,IB+Bdir*i,:] = V[jBA-sign*j,IA+i,:];
            end
        end
    end

    return V
end