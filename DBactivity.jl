#= 
    This uses JuliaDB to manage all of the 
    array data in the program in a non-RAM
    intensive manner. Large datasets are 
    stored in files on disk, and can be
    saved and accessed at any time. 

        !!! This may be written over using
        HDF5. I have never used JuliaDB and
        am just learning about it.
=# 

DBStore = "./DBStore";

function saveDomainData(
    imax, jmax, x, y, xc, yc, xc_g, yc_g, 
    Ai, Aj, ni_x, ni_y, nj_x, nj_y, Volume
)

    dims = [imax, jmax];
    grid = table(x, y, xc, yc, xc_g, yc_g, 
        names = [:x, :y, :xc, :yc, :xc_g, :yc_g]);
    gridinfo = 
        table(dims, Ai, Aj, ni_x, ni_y, nj_x, nj_y, Volume,
        names = [:dims, :Ai, :ni_x, :ni_y, :nj_x, :nj_y, :Volume]);

    savegrid = save(grid, DBStore);
    saveinfo = save(gridinfo, DBStore);

end

function saveInitialData(
    U, S, VL, VR, VB, VT, F, G, Res, MaxSpeed, dt, dt_min, L2
)

    initData = table(U, S, VL, VR, VB, VT, F, G, Res,
        names = [:U, :S, :VL, :VR, :VB, :VT, :F, :G, :Res]);
    initInfo = table(MaxSpeed, dt, dt_min, L2,
        names = [:MaxSpeed, :dt, :dt_min, :L2]);

    savedata = save(initData, DBStore);
    saveinfo = save(initInto, DBStore);

end