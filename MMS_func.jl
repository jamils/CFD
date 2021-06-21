function solveSourceMMS(xc, yc)
    L = 1;
    S = zeros(jmax, imax, 4); # Source terms
    #=
        1 → ρ
        2 → ρu
        3 → ρv
        4 → ρeₜ
    =#

    if supersonic == 0 # subsonic table A.2
        ρ₀ = 1; ρy = 0.15; ρy = -0.1;
        uvel₀ = 70; uvelx = 5; uvely = -7;
        vvel₀ = 90; vvelx = -15; vvely = 8.5;
        wvel₀ = 0; wvelx = 0; wvely = 0;
        p₀ = 1e5; px = 0.2e5; py = 0.5e5;

    elseif supersonic == 1 # supersonic table A.1
        ρ₀ = 1; ρy = 0.15; ρy = -0.1;
        uvel₀ = 800; uvelx = 50; uvely = -30;
        vvel₀ = 800; vvelx = -75; vvely = 40;
        wvel₀ = 0; wvelx = 0; wvely = 0;
        p₀ = 1e5; px = 0.2e5; py = 0.5e5;
    end

    # Apply from mathematica notebook
    for i ∈ 1:imax
        for j ∈ 1:jmax
            x = xc[j,i];
            y = yc[j,i];

            S[j,i,1] = (3*π*uvelx*cos((3*π*x)/(2*L))*(ρ₀ + ρy*cos((π*y)/(2*L)) + 
                ρx*sin((π*x)/L)))/(2*L) + (2*π*vvely*cos((2*π*y)/(3*L))*
                (ρ₀ + ρy*cos((π*y)/(2*L)) + ρx*sin((π*x)/L)))/(3*L) + 
                (π*ρx*cos((π*x)/L)*(uvel₀ + uvely*cos((3*π*y)/(5*L)) + 
                uvelx*sin((3*π*x)/(2*L))))/L - (π*ρy*sin((π*y)/(2*L))*
                (vvel₀ + vvelx*cos((π*x)/(2*L)) + vvely*sin((2*π*y)/(3*L))))/(2*L);

            S[j,i,2] = (3*π*uvelx*cos((3*π*x)/(2*L))*(ρ₀ + ρy*cos((π*y)/(2*L)) + 
                ρx*sin((π*x)/L))*(uvel₀ + uvely*cos((3*π*y)/(5*L)) + 
                uvelx*sin((3*π*x)/(2*L))))/L + (2*π*vvely*cos((2*π*y)/(3*L))*
                (ρ₀ + ρy*cos((π*y)/(2*L)) + ρx*sin((π*x)/L))*
                (uvel₀ + uvely*cos((3*π*y)/(5*L)) + 
                uvelx*sin((3*π*x)/(2*L))))/(3*L) + 
                (π*ρx*cos((π*x)/L)*((uvel₀ + uvely*cos((3*π*y)/(5*L)) + 
                uvelx*sin((3*π*x)/(2*L)))^2))/L - (2*π*px*sin((2*π*x)/L))/L - 
                (π*ρy*(uvel₀ + uvely*cos((3*π*y)/(5*L)) + 
                uvelx*sin((3*π*x)/(2*L)))*sin((π*y)/(2*L))*
                (vvel₀ + vvelx*cos((π*x)/(2*L)) + 
                vvely*sin((2*π*y)/(3*L))))/ (2*L) - 
                (3*π*uvely*(ρ₀ + ρy*cos((π*y)/(2*L)) + 
                ρx*sin((π*x)/L))*sin((3*π*y)/(5*L))*
                (vvel₀ + vvelx*cos((π*x)/(2*L)) + vvely*sin((2*π*y)/(3*L))))/(5*L);

            S[j,i,3] = (π*py*cos((π*y)/L))/L - 
                (π*vvelx*sin((π*x)/(2*L))*(ρ₀ + ρy*cos((π*y)/(2*L)) + 
                ρx*sin((π*x)/L))*(uvel₀ + uvely*cos((3*π*y)/(5*L)) + 
                uvelx*sin((3*π*x)/(2*L))))/(2*L) + 
                (3*π*uvelx*cos((3*π*x)/(2*L))*(ρ₀ + ρy*cos((π*y)/(2*L)) + 
                ρx*sin((π*x)/L))*(vvel₀ + vvelx*cos((π*x)/(2*L)) + 
                vvely*sin((2*π*y)/(3*L))))/(2*L) + 
                (4*π*vvely*cos((2*π*y)/(3*L))*(ρ₀ + ρy*cos((π*y)/(2*L)) + 
                ρx*sin((π*x)/L))*(vvel₀ + vvelx*cos((π*x)/(2*L)) + 
                vvely*sin((2*π*y)/(3*L))))/(3*L) + 
                (π*ρx*cos((π*x)/L)*(uvel₀ + uvely*cos((3*π*y)/(5*L)) + 
                uvelx*sin((3*π*x)/(2*L)))*(vvel₀ + vvelx*cos((π*x)/(2*L)) + 
                vvely*sin((2*π*y)/(3*L))))/L - 
                (π*ρy*sin((π*y)/(2*L))*((vvel₀ + vvelx*cos((π*x)/(2*L)) + 
                vvely*sin((2*π*y)/(3*L)))^2))/(2*L); 

            S[j,i,4] = (uvel₀ + uvely*cos((3*π*y)/(5*L)) + 
                uvelx*sin((3*π*x)/(2*L)))*((-2*π*px*sin((2*π*x)/L))/L + 
                (ρ₀ + ρy*cos((π*y)/(2*L)) + 
                ρx*sin((π*x)/L))*((-2*π*px*sin((2* π*x)/L))/((-1 + gamma)*L*
                (ρ₀ + ρy*cos((π*y)/(2*L)) + ρx*sin((π*x)/L))) + 
                ((3*π*uvelx*cos((3*π*x)/(2*L))*(uvel₀ + uvely*cos((3*π*y)/(5*L)) + 
                uvelx*sin((3*π*x)/(2*L))))/L - 
                (π*vvelx*sin((π*x)/(2*L))*(vvel₀ + vvelx*cos((π*x)/(2*L)) + 
                vvely*sin((2*π*y)/(3*L))))/L)/2 - 
                (π*ρx*cos((π*x)/L)*(p0 + px*cos((2*π*x)/L) + 
                py*sin((π*y)/L)))/((-1 + gamma)*L*((ρ₀ + ρy*cos((π*y)/(2*L)) + 
                ρx*sin((π*x)/L))^2))) + (π*ρx*cos((π*x)/L)*(((wvel₀^2) + 
                ((uvel₀ + uvely*cos((3*π*y)/(5*L)) + uvelx*sin((3*π*x)/(2*L)))^2) + 
                ((vvel₀ + vvelx*cos((π*x)/(2*L)) + vvely*sin((2*π*y)/(3*L)))^2))/2 + 
                (p0 + px*cos((2*π*x)/L) + py*sin((π*y)/L))/((-1 + gamma)*
                (ρ₀ + ρy*cos((π*y)/(2*L)) + ρx*sin((π*x)/L)))))/L) + 
                (3*π*uvelx*cos((3*π*x)/(2*L))*(p0 + px*cos((2*π*x)/L) + 
                py*sin((π*y)/L) + (ρ₀ + ρy*cos((π*y)/(2*L)) + 
                ρx*sin((π*x)/L))*(((wvel₀^2) + ((uvel₀ + uvely*cos((3*π*y)/(5*L)) + 
                uvelx*sin((3*π*x)/(2*L))^2) + ((vvel₀ + vvelx*cos((π*x)/(2*L)) + 
                vvely*sin((2*π*y)/(3*L)))^2))/2 + (p0 + px*cos((2*π*x)/L) + 
                py*sin((π*y)/L))/((-1 + gamma)*(ρ₀ + ρy*cos((π*y)/(2*L)) + 
                ρx*sin((π*x)/L)))))))/(2*L) + (2*π*vvely*cos((2*π*y)/(3*L))*
                (p0 + px*cos((2*π*x)/L) + py*sin((π*y)/L) + 
                (ρ₀ + ρy*cos((π*y)/(2*L)) + ρx*sin((π*x)/L))*(((wvel₀^2) + 
                ((uvel₀ + uvely*cos((3*π*y)/(5*L)) + uvelx*sin((3*π*x)/(2*L)))^2) + 
                ((vvel₀ + vvelx*cos((π*x)/(2*L)) + vvely*sin((2*π*y)/(3*L)))^2))/2 + 
                (p0 + px*cos((2*π*x)/L) + py*sin((π*y)/L))/((-1 + gamma)*
                (ρ₀ + ρy*cos((π*y)/(2*L)) + ρx*sin((π*x)/L))))))/(3*L) + 
                (vvel₀ + vvelx*cos((π*x)/(2*L)) + vvely*sin((2*π*y)/(3*L)))*
                ((π*py*cos((π*y)/L))/L - (π*ρy*sin((π*y)/(2*L))*(((wvel₀^2) + 
                ((uvel₀ + uvely*cos((3*π*y)/(5*L)) + uvelx*sin((3*π*x)/(2*L)))^2) + 
                ((vvel₀ + vvelx*cos((π*x)/(2*L)) + vvely*sin((2*π*y)/(3*L)))^2))/2 + 
                (p0 + px*cos((2*π*x)/L) + py*sin((π*y)/L))/((-1 + gamma)*
                (ρ₀ + ρy*cos((π*y)/(2*L)) + ρx*sin((π*x)/L)))))/(2*L) + 
                (ρ₀ + ρy*cos((π*y)/(2*L)) + ρx*sin((π*x)/L))*
                ((π*py*cos((π*y)/L))/((-1 + gamma)*L*(ρ₀ + ρy*cos((π*y)/(2*L)) + 
                ρx*sin((π*x)/L))) + ((-6*π*uvely*(uvel₀ + uvely*cos((3*π*y)/(5*L)) + 
                uvelx*sin((3*π*x)/(2*L)))*sin((3*π*y)/(5*L)))/(5*L) + 
                (4*π*vvely*cos((2*π*y)/(3*L))*(vvel₀ + vvelx*cos((π*x)/(2*L)) + 
                vvely*sin((2*π*y)/(3*L))))/(3*L))/2 + (π*ρy*sin((π*y)/(2*L))*
                (p0 + px*cos((2*π*x)/L) + py*sin((π*y)/L)))/(2*(-1 + gamma)*L*
                ((ρ₀ + ρy*cos((π*y)/(2*L)) + ρx*sin((π*x)/L))^2))));
        end
    end

    return S
end

function solveSolutionMMS(xc_g, yc_g)
    L = 1;
    V_MMS = zeros(jmax, imax, 4);
    #=
        1 → ρ
        2 → u
        3 → v
        4 → p
    =#

    if supersonic == 0 # subsonic table A.2
        ρ₀ = 1; ρy = 0.15; ρy = -0.1;
        uvel₀ = 70; uvelx = 5; uvely = -7;
        vvel₀ = 90; vvelx = -15; vvely = 8.5;
        wvel₀ = 0; wvelx = 0; wvely = 0;
        p₀ = 1e5; px = 0.2e5; py = 0.5e5;

    elseif supersonic == 1 # supersonic table A.1
        ρ₀ = 1; ρy = 0.15; ρy = -0.1;
        uvel₀ = 800; uvelx = 50; uvely = -30;
        vvel₀ = 800; vvelx = -75; vvely = 40;
        wvel₀ = 0; wvelx = 0; wvely = 0;
        p₀ = 1e5; px = 0.2e5; py = 0.5e5;
    end

    for i ∈ 1:imax
        for j ∈ 1:jmax
            x = xc_g[j,i];
            y = yc_g[j,i];

            V_MMS[j,i,1] = ρ₀ + ρy * cos((π*y)/(2*L)) + ρx * sin((π*x)/L);
            V_MMS[j,i,2] = uvel₀ + uvely * cos((3*π*y)/(5*L)) + uvelx * 
                sin((3*π*x)/(2*L));
            V_MMS[j,i,3] = vvel₀ + vvelx * cos((3*π*y)/(5*L)) + uvelx * 
                sin((3*π*x)/(2*L));
            V_MMS[j,i,4] = p₀ + px * cos((2*π*x)/L) + py * sin((π*y)/L);
        end
    end

    return V_MMS
end

function setBCMMS(V_MMS)
    V = zeros(jmax, imax, 4);

    ni_g = imax;
    nj_g = jmax;
    ni = imax - 2*num_ghost;
    nj = jmax - 2*num_ghost;

    bot = num_ghost - 1;
    top = num_ghost + nj;
    left = num_ghost - 1;
    right = num_ghost + ni;

    for i ∈ num_ghost:(ni+num_ghost)
        for j ∈ 1:num_ghost
            @. V[bot-j,i,:] = V_MMS[bot-j,i,:];
            @. V[top+j,i,:] = V_MMS[top+j,i,:];
        end
    end
       
    for i ∈ 1:num_ghost
        for j ∈ num_ghost:(num_ghost+nj)
            @. V[j,left-i,:] = V_MMS[j,left-i,:];
            @. V[j,right+i,:] = V_MMS[j,right+i,:];
        end
    end

    return V
end