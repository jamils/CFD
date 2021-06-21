function initialize(xc_g, yc_g)
    if case == 1 # Curvilinear
        if supersonic == 0 # Subsonic
            ρ = 1; # [kg/m³]
            u = 70; # [m/s]
            v = 90; # [m/s]
            p = 1e5; # [Pa]

        elseif supersonic == 1 # Supersonic
            ρ = 1; # [kg/m³]
            u = 800; # [m/s]
            v = 800; # [m/s]
            p = 1e5; # [Pa]
        end

    elseif case == 2 # 30ᵒ Inlet
        M = 4.0;
        T = 217; # [K]
        p = 12270; # [Pa]
        ρ = p/(R*T); # [kg/m³]
        u = M*√(γ*R*T); # [m/s]
        v = 0;

    elseif case == 3 # NACA Airfoil
        if AngleOfAttack == 1 # 0ᵒ
            M = 0.84;
            p = 65855.8; # [Pa]
            AOA = 0;

        elseif AngleOfAttack == 2 # 8ᵒ
            M = 0.75;
            p = 67243.5; # [Pa]
            AOA = 8;
        end

        T = 300; # [K]
        ρ = p/(R*T);
        u = M*√(γ*R*T) * cos(AOA*π/180);
        v = M*√(γ*R*T) * sin(AOA*π/180);

    end

    # Create Primitive variables
    V = zeros(jmax, imax, 4);
    @. V[:,:,1] = ρ;
    @. V[:,:,2] = u;
    @. V[:,:,3] = v;
    @. V[:,:,4] = p;

    return V
end