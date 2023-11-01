
"""
    get_gm(id)

Return the gravitational parameter of a body with the given `id`.
"""
function get_gm(id)
    if id == 10
        return 132712440017.99
    elseif id == 301
        return 4902.8005821478
    elseif id == 399
        return 398600.4415
    else
        error("Unknown target id: $id")
    end
end