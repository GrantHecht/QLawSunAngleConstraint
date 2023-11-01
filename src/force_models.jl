
"""
    ForceModel

Struct represnting the force model used by the `DSMManeuverPlanner.jl` package.

# Fields
- `epoch0::Float64`: Initial epoch in SPICE ET (TBD)
- `mu::Float64`: Primary bodies gravitational parameter
- `perterbations::Bool`: If false, no perturbations are included in force model
- `third_body::Bool`: If true, third body perturbations are included in force model
- `ephem::ET`: Ephemerides for third bodies
"""
struct ForceModel{ET <: Union{Nothing, Ephemerides}}
    # Initial epoch in SPICE ET (TBD)
    epoch0::Float64

    # Primary bodies gravitational parameter
    mu::Float64

    # Perturbation flags
    perterbations::Bool     # If false, no perturbations are included in force model
    third_body::Bool        # If true, third body perturbations are included in force model

    # Ephemerides for third bodies
    ephem::ET
end

"""
    ForceModel(epoch0, cb_id; perterbations = true, third_body = true)

Construct a `ForceModel` struct for the given central body.

# Arguments
- `epoch0::AbstractFloat`: Initial epoch in SPICE ET (TBD)
- `cb_id::Int`: NAIF ID of central body

# Keywords
- `perterbations::Bool=true`: If false, no perturbations are included in force model
- `third_body::Bool=true`: If true, third body perturbations are included in force model

# Returns
- `p::ForceModel`: The constructed `ForceModel` struct

# Notes
- Currently may call SPICE. This should be addressed in the future.
- Currently hardcoding μ for Enceladus because the value provided to GMAT is different than
    what is provided in gm_de440.tpc
- The third body perturbations will include the gravitational effects of each target
    body stored in the serialized ephemrides file. This means that if you want to include
    the gravitational effects of a body that is not a target, you must add it to the
    serialized ephemerides file when you create it.
"""
function ForceModel(
    epoch0::AbstractFloat, cb_id::Int;
    perterbations = true, third_body = true,
)
    # Furnsh kernels
    furnsh_kernals()

    # Get central bodies gravitational parameter
    mu = get_gm(cb_id)

    # Unload kernels
    unload_kernals()

    # Get third body ephemerides if desired
    ephem = third_body ? load_serialized_ephemerides() : nothing

    return ForceModel(epoch0, mu, perterbations, third_body, ephem)
end

"""
    ForceModel(epoch0, mu; perterbations = true, third_body = true)

Construct a `ForceModel` struct for the given central body.

# Arguments
- `epoch0::AbstractFloat`: Initial epoch in SPICE ET (TBD)
- `mu::AbstractFloat`: Primary bodies gravitational parameter

# Keywords
- `perterbations::Bool=true`: If false, no perturbations are included in force model
- `third_body::Bool=true`: If true, third body perturbations are included in force model

# Returns
- `p::ForceModel`: The constructed `ForceModel` struct

# Notes
- The third body perturbations will include the gravitational effects of each target
    body stored in the serialized ephemrides file. This means that if you want to include
    the gravitational effects of a body that is not a target, you must add it to the
    serialized ephemerides file when you create it.
"""
function ForceModel(
    epoch0::AbstractFloat, mu::AbstractFloat;
    perterbations = true, third_body = true,
)
    # Get third body ephemerides if desired
    ephem = load_serialized_ephemerides() 

    return ForceModel(epoch0, mu, perterbations, third_body, ephem)
end

"""
    get_perturbing_accelerations(p::ForceModel, state, t)

Compute the perturbing accelerations in the S/Cs RTN frame for the given state and time.

# Arguments
- `p::ForceModel`: The force model struct
- `state::AbstractVector`: The MEE state vector
- `t::AbstractFloat`: The current time in seconds past the initial epoch.

# Returns
- `ap::SVector{3,Float64}`: The perturbing accelerations in the S/Cs RTN frame
"""
function get_perturbing_accelerations(p::ForceModel, state, t)
    # Return null acceleration if we're not including perts
    if p.perterbations == false
        return SVector(0.0, 0.0, 0.0)
    end

    # Accumulate perturbing accelerations
    ar = 0.0
    at = 0.0
    an = 0.0
    if p.third_body == true
        tbp = third_body_perturbations(p, state, t)
        ar += tbp[1]
        at += tbp[2]
        an += tbp[3]
    end

    # Return accelerations
    return SVector(ar, at, an)
end

"""
    get_perturbing_accelerations_and_partials(p::ForceModel, state, t)

Compute the perturbing accelerations and partial derivatives of the perturbing accelerations
    w.r.t the MEE state.

# Arguments
- `p::ForceModel`: The force model struct
- `state::AbstractVector`: The MEE state vector
- `t::AbstractFloat`: The current time in seconds past the initial epoch.

# Returns
- `ap::SVector{3,Float64}`: The perturbing accelerations in the S/Cs RTN frame
- `dapdx::SMatrix{3,6,Float64}`: The partial derivatives of the perturbing accelerations
    w.r.t the MEE state
"""
function get_perturbing_accelerations_and_partials(p::ForceModel, state, t)
    # Return null acceleration and partials if we're not including perts
    if p.perterbations == false
        ap      = SVector(0.0, 0.0, 0.0)
        dapdx   = SMatrix{3,6}(
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        )
        return ap, dapdx
    end

    # Accumulate perturbing accelerations and partials
    ar  = 0.0
    at  = 0.0
    an  = 0.0
    # Note: Choosing to use MutableArray here for accumulating the partials
    # but more performant options may exist
    dapdx_mutable = MMatrix{3,6}(
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    )
    if p.third_body == true
        tbp, dtbpdx = third_body_perturbations_and_partials(p, state, t)
        ar += tbp[1]
        at += tbp[2]
        an += tbp[3]
        dapdx_mutable += dtbpdx
    end

    # Return accelerations and partials
    ap = SVector(ar, at, an)
    dapdx = SMatrix{3, 6}(dapdx_mutable)
    return ap, dapdx
end

"""
    third_body_perturbations(p::ForceModel, state, t)

Compute the third-body perturbing accelerations in the S/Cs RTN frame for the given state
    and time.

# Arguments
- `p::ForceModel`: The force model struct
- `state::AbstractVector`: The MEE state vector
- `t::AbstractFloat`: The current time in seconds past the initial epoch.

# Returns
- `ap::SVector{3,Float64}`: The perturbing accelerations in the S/Cs RTN frame
"""
function third_body_perturbations(p::ForceModel{Ephemerides}, mee_state, t)
    # Convert state to cartesian
    cart = convert_state(mee_state, MEE, Cartesian, p.mu)
    pos  = SVector(cart[1], cart[2], cart[3])

    # Shift sim time to SPICE ET
    t_et = p.epoch0 + t

    # Get gravitational parameters
    GMs = get_gravitational_parameter(p.ephem)

    # Get third body positions at time t_et
    tbpos = interpolate(p.ephem, t_et)

    # Compute accelerations in inertial frame
    ax = 0.0
    ay = 0.0
    az = 0.0
    @inbounds for i in eachindex(tbpos)
        # Get position of third body
        pos_tb  = SVector(tbpos[i][1], tbpos[i][2], tbpos[i][3])

        # Compute acceleration due to ith body
        #a_tb    = _third_body(pos, pos_tb, GMs[i])
        a_tb    = _third_body_battin(pos, pos_tb, GMs[i])

        # Accumulate accelerations
        ax     += a_tb[1]
        ay     += a_tb[2]
        az     += a_tb[3]
    end

    # Rotate into RTN frame
    a_inertial  = SVector(ax, ay, az)
    a_rtn       = _inertial_to_rtn_vec(a_inertial, cart)

    # Return perturbation in RTN frame
    return a_rtn
end
function third_body_perturbations(p::ForceModel{Nothing}, state, t)
    return SVector(0.0, 0.0, 0.0)
end

"""
    third_body_perturbations_and_partials(p::ForceModel, state, t)

Compute the third-body perturbing accelerations and partial derivatives of the
    third-body perturbing accelerations w.r.t the MEE state.

# Arguments
- `p::ForceModel`: The force model struct
- `state::AbstractVector`: The MEE state vector
- `t::AbstractFloat`: The current time in seconds past the initial epoch.

# Returns
- `ap::SVector{3,Float64}`: The perturbing accelerations in the S/Cs RTN frame
- `dapdx::SMatrix{3,6,Float64}`: The partial derivatives of the perturbing accelerations
    w.r.t the MEE state
"""
function third_body_perturbations_and_partials(p::ForceModel{Ephemerides}, mee_state, t)
    # Convert state to cartesian
    cart = convert_state(mee_state, MEE, Cartesian, p.mu)
    pos  = SVector(cart[1], cart[2], cart[3])

    # Compute state conversion jacobian
    dcartdmee = convert_state_partials(mee_state, MEE, Cartesian, p.mu)

    # Shift sim time to SPICE ET
    t_et = p.epoch0 + t

    # Get gravitational parameters
    GMs = get_gravitational_parameter(p.ephem)

    # Get third body positions at time t_et
    tbpos = interpolate(p.ephem, t_et)

    # Compute accelerations in inertial frame
    ax = 0.0
    ay = 0.0
    az = 0.0
    dadcart = MMatrix{3,6}(
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    )
    @inbounds for i in eachindex(tbpos)
        # Get position of third body
        pos_tb  = SVector(tbpos[i][1], tbpos[i][2], tbpos[i][3])

        # Compute acceleration due to ith body
        a_tb, a_tb_partials = _third_body_battin_and_partials(pos, pos_tb, GMs[i])

        # Accumulate accelerations
        ax     += a_tb[1]
        ay     += a_tb[2]
        az     += a_tb[3]

        # Accumulate partials
        dadcart += a_tb_partials
    end

    # Rotate into RTN frame
    a_inertial  = SVector(ax, ay, az)
    a_rtn       = _inertial_to_rtn_vec(a_inertial, cart)

    # Construct final partial
    dartndcart  = _inertial_to_rtn_vec_partial_cart(a_inertial, dadcart, cart)
    dartndmee   = dartndcart*dcartdmee

    # Return perturbation in RTN frame
    return a_rtn, dartndmee
end
function third_body_perturbations_and_partials(p::ForceModel{Nothing}, state, t)
    ap = SVector(0.0, 0.0, 0.0)
    dapdx = SMatrix{3,6}(
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    )
    return ap, dapdx
end

"""
    _third_body(pos, pos_tb, gm_tb)

Compute the third body acceleration in the inertial frame using the standard formulation
    for accelerations due to third bodies

# Arguments
- `pos::SVector{3,Float64}`: Position of spacecraft in inertial frame
- `pos_tb::SVector{3,Float64}`: Position of third body in inertial frame
- `gm_tb::Float64`: Gravitational parameter of third body

# Returns
- `a_tb::SVector{3,Float64}`: Third body acceleration in inertial frame
"""
function _third_body(pos, pos_tb, gm_tb)
    # Compute vector from sc to third body
    pos_tb_sc = SVector(
        pos_tb[1] - pos[1],
        pos_tb[2] - pos[2],
        pos_tb[3] - pos[3],
    )

    # Compute inverse norm cubed of both position vectors
    inv_norm_pos_tb_sc  = 1.0 / norm(pos_tb_sc)
    inv_norm_pos_tb     = 1.0 / norm(pos_tb)
    inv_norm_pos_tb_sc3 = inv_norm_pos_tb_sc*inv_norm_pos_tb_sc*inv_norm_pos_tb_sc
    inv_norm_pos_tb3    = inv_norm_pos_tb*inv_norm_pos_tb*inv_norm_pos_tb

    # Compute third body acceleration
    return SVector(
        gm_tb*(pos_tb_sc[1]*inv_norm_pos_tb_sc3 - pos_tb[1]*inv_norm_pos_tb3),
        gm_tb*(pos_tb_sc[2]*inv_norm_pos_tb_sc3 - pos_tb[2]*inv_norm_pos_tb3),
        gm_tb*(pos_tb_sc[3]*inv_norm_pos_tb_sc3 - pos_tb[3]*inv_norm_pos_tb3),
    )
end

"""
    _third_body_battin(pos, pos_tb, gm_tb)

Compute the third body acceleration in the inertial frame using Battin's method

# Arguments
- `pos::SVector{3,Float64}`: Position of spacecraft in inertial frame
- `pos_tb::SVector{3,Float64}`: Position of third body in inertial frame
- `gm_tb::Float64`: Gravitational parameter of third body

# Returns
- `a_tb::SVector{3,Float64}`: Third body acceleration in inertial frame
"""
function _third_body_battin(pos, pos_tb, gm_tb)
    # Compute relative position of spacecraft with respect to the third body
    d   = SVector(pos[1] - pos_tb[1], pos[2] - pos_tb[2], pos[3] - pos_tb[3])

    # Compute Battin's F(q) function
    Fqk = _compute_fk(pos, pos_tb)

    # Compute acceleration
    dk3     = sqrt(d[1]*d[1] + d[2]*d[2] + d[3]*d[3])^3
    dk3_inv = 1.0 / dk3
    return SVector(
        -gm_tb*dk3_inv*(pos[1] + Fqk*pos_tb[1]),
        -gm_tb*dk3_inv*(pos[2] + Fqk*pos_tb[2]),
        -gm_tb*dk3_inv*(pos[3] + Fqk*pos_tb[3]),
    )
end

"""
    _third_body_battin_and_partials(pos, pos_tb, gm_tb)

Compute the third body acceleration and its partial derivative w.r.t. the cartisian state.

# Arguments
- `pos::SVector{3,Float64}`: Position of spacecraft in inertial frame
- `pos_tb::SVector{3,Float64}`: Position of third body in inertial frame
- `gm_tb::Float64`: Gravitational parameter of third body

# Returns
- `atb::SVector{3,Float64}`: Third body acceleration in inertial frame
- `datbdx::SMatrix{3,6,Float64}`: Partial derivative of third body acceleration
    w.r.t. the cartisian state
"""
function _third_body_battin_and_partials(pos, pos_tb, gm_tb)
    # Compute relative position of spacecraft with respect to the third body
    d   = SVector(pos[1] - pos_tb[1], pos[2] - pos_tb[2], pos[3] - pos_tb[3])

    # Compute Battin's F(q) function
    Fqk = _compute_fk(pos, pos_tb)

    # Compute acceleration
    dk      = sqrt(d[1]*d[1] + d[2]*d[2] + d[3]*d[3])
    dk3     = dk^3
    dk3_inv = 1.0 / dk3
    rpfs    = SVector(
        pos[1] + Fqk*pos_tb[1],
        pos[2] + Fqk*pos_tb[2],
        pos[3] + Fqk*pos_tb[3],
    )
    atb     = -gm_tb*dk3_inv*rpfs

    # Compute partials
    dk5     = dk3*dk*dk
    dk5_inv = 1.0 / dk5
    term1   = 3.0*gm_tb*dk5_inv*rpfs*transpose(SVector(d[1], d[2], d[3], 0.0, 0.0, 0.0))

    drdx    = SMatrix{3,6}(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
    )
    term2   = -gm_tb*dk3_inv*(drdx + pos_tb*_compute_fk_partial_x(pos, pos_tb))

    datbdx  = term1 + term2
    return atb, datbdx
end


"""
    _compute_qk(r, s)

Compute the qk parameter for Battin's method

# Arguments
- `r::AbstractVector`: Position of spacecraft in inertial frame
- `s::AbstractVector`: Position of third body in inertial frame

# Returns
- `qk::Float64`: qk parameter
"""
function _compute_qk(r, s)
    sds = @inbounds s[1]*s[1] + s[2]*s[2] + s[3]*s[3]
    qk  = @inbounds (r[1]*(r[1] - 2.0*s[1]) +
                     r[2]*(r[2] - 2.0*s[2]) +
                     r[3]*(r[3] - 2.0*s[3])) / sds
    return qk
end

"""
    _compute_qk_partial_x(r, s)

Compute the partial derivative of qk w.r.t. the cartisian state.

# Arguments
- `r::AbstractVector`: Position of spacecraft in inertial frame
- `s::AbstractVector`: Position of third body in inertial frame

# Returns
- `dqkdx::SVector{6,Float64}`: Partial derivative of qk w.r.t. the cartisian state.

# Notes
This function is partially generated by `src/code_gen/q_partials.m`.
"""
function _compute_qk_partial_x(r, s)
    t2 = s[1]*s[1]
    t3 = s[2]*s[2]
    t4 = s[3]*s[3]
    t5 = t2+t3+t4
    t6 = 1.0/t5
    t7 = 2.0*t6
    return transpose(SVector(
        t7*(r[1]-s[1]),
        t7*(r[2]-s[2]),
        t7*(r[3]-s[3]),
        0.0,
        0.0,
        0.0,
    ))
end

"""
    _compute_fk(qk)

Compute the fk parameter for Battin's method given qk

# Arguments
- `qk::Float64`: qk parameter

# Returns
- `fk::Float64`: fk parameter
"""
_compute_fk(qk) = qk*(3.0 + 3.0*qk + qk*qk) / (1.0 + sqrt(1.0 + qk)^3)

"""
    _compute_fk(r, s)

Compute the fk parameter for Battin's method

# Arguments
- `r::AbstractVector`: Position of spacecraft in inertial frame
- `s::AbstractVector`: Position of third body in inertial frame

# Returns
- `fk::Float64`: fk parameter
"""
_compute_fk(r, s) = _compute_fk(_compute_qk(r, s))

"""
    _compute_fk_partial_qk(qk)

Compute the partial derivative of fk w.r.t. qk

# Arguments
- `qk::Float64`: qk parameter

# Returns
- `dfdq::Float64`: Partial derivative of fk w.r.t. qk

# Notes
This function is partially generated by `src/code_gen/Fq_partials.m`.
"""
function _compute_fk_partial_qk(qk)
    t2 = qk*3.0;
    t3 = qk+1.0;
    t4 = qk*qk;
    t5 = t3^(3.0/2.0);
    t6 = t2+t4+3.0;
    t7 = t5+1.0;
    t8 = 1.0/t7;
    return t6*t8+qk*t8*(t3*2.0+1.0)-qk*sqrt(t3)*t6*(t8*t8)*(3.0/2.0);
end

"""
    _compute_fk_partial_x(r, s)

Compute the partial derivative of fk w.r.t. the cartisian state.

# Arguments
- `r::AbstractVector`: Position of spacecraft in inertial frame
- `s::AbstractVector`: Position of third body in inertial frame

# Returns
- `dfdx::SVector{6,Float64}`: Partial derivative of fk w.r.t. the cartisian state.
"""
function _compute_fk_partial_x(r, s)
    dfdq = _compute_fk_partial_qk(_compute_qk(r, s))
    dqdx = _compute_qk_partial_x(r, s)
    return dfdq*dqdx
end
