function compute2DFlux(ni_x, ni_y, nj_x, nj_y, VL, VR, VB, VT)
    if upwind == 1
        F = computeFluxVL(ni_x, ni_y, VL, VR);
        G = computeFluxVL(nj_x, nj_y, VB, VT);

    elseif upwind == 2
        F = computeFluxRoe(ni_x, ni_y, VL, VR);
        G = computeFluxRoe(nj_x, nj_y, VB, VT);
    end

    return F, G
end

function computeFluxVL(nx, ny, VL, VR)
    Fc = zeros(4);
    Fp = zeros(4);
    F = zeros(jmax, imax, 4);

    for i ∈ 1:imax
        for j ∈ 1:jmax
            # Sound speed
            aL = √(γ*VL[j,i,4] / VL[j,i,1]);
            aR = √(γ*VR[j,i,4] / VR[j,i,1]);

            println("aL = $aL")
            println("aR = $aR")

            # Velocity
            uL = VL[j,i,2] * nx[j,i] + VL[j,i,3] * ny[j,i];
            uR = VR[j,i,2] * nx[j,i] + VR[j,i,3] * ny[j,i];

            # Mach number
            ML = uL / aL;
            MR = uR / aR;

            println("nx[$j,$i] = $(nx[j,i])")
            println("ny[$j,$i] = $(ny[j,i])")

            M⁺ = 0.25 * (ML + 1)^2;
            M⁻ = 0.25 * (MR - 1)^2;

            βL = - max(0, 1 - Int64(round(abs(ML))));
            βR = - max(0, 1 - Int64(round(abs(MR))));

            α⁺ = 0.5 * (1 + copysign(1, ML));
            α⁻ = 0.5 * (1 - copysign(1, MR));

            c⁺ = α⁺ * (1 + βL) * ML - βL * M⁺;
            c⁻ = α⁻ * (1 + βR) * MR - βR * M⁻;

            # Average pressure
            p̄⁺ = M⁺ * (-ML + 2);
            p̄⁻ = M⁻ * (-MR - 2);

            D⁺ = α⁺ * (1 + βL) - βL * p̄⁺;
            D⁻ = α⁻ * (1 + βR) - ΒR * p̄⁻;

            # Enthalpy
            hₜL = (γ/(γ-1)) * VL[j,i,4] / VL[j,i,1] + 
                0.5 * (VL[j,i,2]^2 + VL[j,i,3]^2);
            hₜR = (γ/(γ-1)) * VR[j,i,4] / VR[j,i,1] + 
                0.5 * (VR[j,i,2]^2 + VR[j,i,3]^2);
            
            # Compute convect flux
            Fc[1] = VL[j,i,1] * aL * c⁺ * VR[j,i,1] * aR * c⁻;

            Fc[2] = VL[j,i,1] * aL * c⁺ * VL[j,i,2] + 
                VR[j,i,1] * aR * c⁻ * VR[j,i,2];

            Fc[3] = VL[j,i,1] * aL * c⁺ * VL[j,i,3] + 
                VR[j,i,1] * aR * c⁻ * VR[j,i,3];

            Fc[4] = VL[j,i,1] * aL * c⁺ * hₜL + VR[j,i,1] * aR * c⁻ * hₜR;

            # Compute pressure flux
            Fp[1] = 0;

            Fp[2] = D⁺ * nx[j,i] * VL[j,i,4] + D⁻ * nx[j,i] * VR[j,i,4];
            
            Fp[3] = D⁺ * ny[j,i] * VL[j,i,4] + D⁻ * ny[j,i] * VR[j,i,4];

            Fp[4] = 0;

            @. F[j,i,:] = Fc + Fp;

        end
    end

    return F
end

function computeFLuxRoe(nx, ny, VL, VR)
    Ε = 0.1; # correction

    r1 = zeros(4);
    r2 = zeros(4);
    r3 = zeros(4);
    r4 = zeros(4);

    FL = zeros(4);
    FR = zeros(4);

    for j ∈ 1:jmax
        for i ∈ 1:imax
            R = √(VR[j,i,1] / VL[j,i,1]);
            ρ = R * VL[j,i,1];
            u = (R * VR[j,i,2] + VL[j,i,2]) / (R+1);
            v = (R * VR[j,i,3] + VL[j,i,3]) / (R+1);
            Û = u * nx[j,i] + v*ny[j,i];

            hₜL = (γ/(γ-1)) * VL[j,i,4] / VL[j,i,1] + 
                0.5 * (VL[j,i,2]^2 + VL[j,i,3]^2);
            hₜR = (γ/(γ-1)) * VR[j,i,4] / VR[j,i,1] +
                0.5 * (VR[j,i,2]^2 + VR[j,i,3]^2);
            hₜ = (R*hₜR + hₜL) / (R+1);

            a = √((γ-1) * (hₜ - 0.5 * (u^2 + v^2)));

            λ₁ = abs(Û);
            λ₂ = abs(Û);
            λ₃ = abs(Û + a);
            λ₄ = abs(Û - a);

            if λ₁ < (2*Ε*a)
                λ₁ = λ₁^2 / (4*Ε*a) + Ε*a;
                λ₂ = λ₂^2 / (4*Ε*a) + Ε*a;
            end

            if λ₃ < (2*Ε*a)
                λ₃ = λ₃^2 / (4*Ε*a) + Ε*a;
            end

            if λ₄ < (2*Ε*a)
                λ₄ = λ₄^2 / (4*Ε*a) + Ε*a;
            end

            r1[1] = 1;
            r1[2] = u;
            r1[3] = v;
            r1[4] = 0.5 * (u^2 + v^2);

            r2[1] = 0;
            r2[2] = ρ * ny[j,i];
            r2[3] = -ρ * nx[j,i];
            r2[4] = ρ * (u*ny[j,i] - v*nx[j,i]);

            r3[1] = 0.5 * ρ/a;
            r3[2] = (0.5 * ρ/a) * (u + a*nx[j,i]);
            r3[3] = (0.5 * ρ/a) * (v + a*ny[j,i]);
            r3[4] = (0.5 * ρ/a) * (hₜ + a*Û);

            r4[1] = -0.5 * ρ/a;
            r4[2] = (-0.5 * ρ/a) * (u - a*nx[j,i]);
            r4[3] = (-0.5 * ρ/a) * (v - a*ny[j,i]);
            r4[4] = (-0.5 * ρ/a) * (hₜ - a*Û);

            Δρ = VR[j,i,1] - VL[j,i,1];
            Δu = VR[j,i,2] - VL[j,i,2];
            Δv = VR[j,i,3] - VL[j,i,3];
            Δp = VR[j,i,4] - VL[j,i,4];

            # Characteristic variables
            dw1 = Δρ - Δp / (a^2);
            dw2 = Δu * ny[j,i] - Δv * nx[j,i];
            dw3 = Δu * nx[j,i] + Δv * ny[j,i] + Δp / (ρ*a);
            dw4 = Δu * nx[j,i] + Δv * ny[j,i] - Δp / (ρ*a);

            # Convert to 1D
            ÛL = (VL[j,i,2]*nx[j,i] + VL[j,i,3]*ny[j,i]);
            ÛR = (VR[j,i,2]*nx[j,i] + VR[j,i,3]*ny[j,i]);

            FL[1] = VL[j,i,1] * ÛL;
            FR[1] = VR[j,i,1] * ÛR;

            FL[2] = VL[j,i,1] * VL[j,i,2] * ÛL + VL[j,i,4]*nx[j,i];
            FR[2] = VR[j,i,1] * VR[j,i,2] * ÛR + VR[j,i,4]*nx[j,i];

            FL[3] = VL[j,i,1] * VL[j,i,3] * ÛL + VL[j,i,4]*ny[j,i];
            FR[3] = VR[j,i,1] * VR[j,i,3] * ÛR + VR[j,i,4]*ny[j,i];

            FL[4] = VL[j,i,1] * hₜL * ÛL;
            FR[4] = VR[j,i,1] * hₜR * ÛR;

            sum2 = zeros(4);

            # Combine first and second order
            @. sum2 = 0.5 * (abs(λ₂)*dw1*r1 + abs(λ₂)*dw2*r2 +
                abs(λ₃)*dw3*r3 + abs(λ₄)*dw4*r4);

            @. F[j,i,:] = 0.5 * (FL + FR) - sum2;

        end
    end

    return F
end