using Base: Float64
using DelimitedFiles

path = "Project_Files/Grids/Inlet.33x17.grd";

stream = readdlm(path);

ni = stream[2,1];
nj = stream[2,2];
nk = stream[2,3];

xc = zeros(nj, ni);
yc = zeros(nj, ni);
zc = zeros(nj, ni);

xn = zeros(nj+1, ni+1);
yn = zeros(nj+1, ni+1);
zn = zeros(nj+1, ni+1);

global x_index = 4;
global y_index = x_index + (nj+1) * (ni+1) * (nk+1);
global z_index = y_index + (nj+1) * (ni+1) * (nk+1);

xdim, ydim = size(stream);
len = length(stream);
data = Array{Float64}(undef, len);

for i ∈ 1:len
    if stream[i] == ""
        stream[i] = 0;
    end
end

for i ∈ 1:xdim
    data[i:(i+3)] = stream[i,:];
end

for i ∈ 1:nj+1
    for j ∈ 1:ni+1
        global x_index;
        global y_index;
        global z_index;
        xn[i,j] = data[x_index];
        yn[i,j] = data[y_index];
        zn[i,j] = data[z_index];
        x_index = x_index + 1;
        y_index = y_index + 1;
        z_index = z_index + 1;
    end
end

for i ∈ 1:nj+1
    for j ∈ 1:ni+1
        xc[i,j] = 0.25 * (xn[i,j] + xn[i+1,j] + xn[i,j+1] + xn[i+1,j+1]);
        yc[i,j] = 0.25 * (yn[i,j] + yn[i+1,j] + yn[i,j+1] + yn[i+1,j+1]);
        zc[i,j] = 0.25 * (zn[i,j] + zn[i+1,j] + zn[i,j+1] + zn[i+1,j+1]);
    end
end