#= 
Reads in data from the FORTRAN function output files,
    then saves is in a simpler Julia format. 
This is only needed to do once for each geometry.
=#

function fread(dims)
    imax, jmax = dims;

    xin = readdlm("x.txt");
    yin = readdlm("y.txt");

    # global imax = 53; # Columns
    # global jmax = 17; # Rows

    x = zeros(jmax, imax);
    y = zeros(jmax, imax);

    for i ∈ imax
        for j ∈ jmax
            x[j,i] = parse(Float64, readline("x.txt"));
            y[j,i] = parse(Float64, readline("y.txt"));
        end
    end

    return x, y
end

function HDF5Save(x, y, dims)
    imax, jmax = dims;
    fid = h5open("griddata.h5", "cw");
    grouping = "$imax x $jmax";
    create_group(fid, grouping);
    g = fid[grouping];
    g["params"] = [imax, jmax];
    g["x"] = x;
    g["y"] = y;
    close(fid);
end

using HDF5
using DelimitedFiles

dims = [53, 17];
x, y = fread(dims);

HDF5Save(x, y, dims);