function findNaN(ℵ)
    nanloc = zeros(1);
    for i ∈ eachindex(ℵ)
        if isnan(ℵ[i])
            append!(nanloc, i);
        end
    end

    if length(nanloc) == 1 && nanloc[end] == 0
        return println("False")
    else
        return nanloc
    end
end

function findzero(ℵ)
    zeroloc = ones(1);
    for i ∈ eachindex(ℵ)
        if iszero(ℵ[i])
            append!(zeroloc, i);
        end
    end

    if length(zeroloc) == 1 && zeroloc[end] == 1
        return println("False")
    else
        return zeroloc
    end
end