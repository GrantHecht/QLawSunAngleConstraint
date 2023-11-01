
"""
    OrbitalStateRepresentation

An abstract orbital state representation flag.
"""
abstract type OrbitalStateRepresentation end

"""
    Cartesian

A Cartesian orbital state representation flag.
"""
struct Cartesian <: OrbitalStateRepresentation end

"""
    Keplerian

A Keplerian orbital state representation flag.
"""
struct Keplerian <: OrbitalStateRepresentation end

"""
    MEE

A Modified Equinoctual Elements (MEE) orbital state representation flag.
"""
struct MEE <: OrbitalStateRepresentation end

"""
    convert_state(
        state::AbstractArray,
        from_state::Type{OrbitalStateRepresentation},
        to_state::Type{OrbitalStateRepresentation},
        mu::Real,
    )

Converts the `state` from `from_state` to `to_state` using the gravitational parameter `mu`.

# Arguments
- `state::AbstractArray`: the state we want to convert
- `from_state::Type{OrbitalStateRepresentation}`: the orbital state representation we want
    to convert state from.
- `to_state::Type{OrbitalStateRepresentation}`: the orbital state representation we want to
    convert state to
- `mu::Real`: the gravitational parameter for the central body

# Returns
- `SVector{6,AbstractFloat}`: the state converted to to_state

# Throws
- `ArgumentError`: Throws error if convert_state method is not defined for combination
    of from_state and to_state
- `ArgumentError`: Throws error if `state` is not of length 6.
- `DomainError`: Throws error if state is outside of the domain of to_state. This occurs
    if state is singular in the to_state representation (i.e., state is a Cartesian state
    representation for a spacecraft on a circular orbit and to_state is Keplerian)
"""
function convert_state(
    state::AbstractArray,
    from_state::Type{OrbitalStateRepresentation},
    to_state::Type{OrbitalStateRepresentation},
    mu::Real
)
    throw(ArgumentError("State representation conversion not implemented."))
    return nothing
end
function convert_state(
    state::AbstractArray,
    from_state::Type{Cartesian},
    to_state::Type{Keplerian},
    mu::Real
)
    length(state) == 6 || throw(ArgumentError("State must be of length 6."))
    return _cart_to_kep(state, mu)
end
function convert_state(
    state::AbstractArray,
    from_state::Type{Keplerian},
    to_state::Type{Cartesian},
    mu::Real
)
    length(state) == 6 || throw(ArgumentError("State must be of length 6."))
    return _kep_to_cart(state, mu)
end
function convert_state(
    state::AbstractArray,
    from_state::Type{Cartesian},
    to_state::Type{MEE},
    mu::Real
)
    length(state) == 6 || throw(ArgumentError("State must be of length 6."))
    return _cart_to_mee(state, mu)
end
function convert_state(
    state::AbstractArray,
    from_state::Type{MEE},
    to_state::Type{Cartesian},
    mu::Real
)
    length(state) == 6 || throw(ArgumentError("State must be of length 6."))
    return _mee_to_cart(state, mu)
end
function convert_state(
    state::AbstractArray,
    from_state::Type{Keplerian},
    to_state::Type{MEE},
    mu::Real
)
    length(state) == 6 || throw(ArgumentError("State must be of length 6."))
    return _kep_to_mee(state, mu)
end
function convert_state(
    state::AbstractArray,
    from_state::Type{MEE},
    to_state::Type{Keplerian},
    mu::Real
)
    length(state) == 6 || throw(ArgumentError("State must be of length 6."))
    return _mee_to_kep(state, mu)
end

"""
    convert_state_partials(
        state::AbstractArray,
        from_state::Type{OrbitalStateRepresentation},
        to_state::Type{OrbitalStateRepresentation},
        mu::Real,
    )

Computates partial derivative of conversion of the `state` from `from_state` to `to_state`
using the gravitational parameter `mu` with respect to `state`.

# Arguments
- `state::AbstractArray`: the state we want to convert
- `from_state::Type{OrbitalStateRepresentation}`: the orbital state representation we want
    to convert state from.
- `to_state::Type{OrbitalStateRepresentation}`: the orbital state representation we want to
    convert state to
- `mu::Real`: the gravitational parameter for the central body

# Returns
- `SMatrix{6,6,AbstractFloat}`: the jacobian of the state conversion to to_state

# Throws
- `ArgumentError`: Throws error if convert_state method is not defined for combination
    of from_state and to_state
- `ArgumentError`: Throws error if `state` is not of length 6.
- `DomainError`: Throws error if state is outside of the domain of to_state. This occurs
    if state is singular in the to_state representation (i.e., state is a Cartesian state
    representation for a spacecraft on a circular orbit and to_state is Keplerian)
"""
function convert_state_partials(
    state::AbstractArray,
    from_state::Type{OrbitalStateRepresentation},
    to_state::Type{OrbitalStateRepresentation},
    mu::Real,
)
    throw(ArgumentError("State representation conversion not implemented."))
    return nothing
end
function convert_state_partials(
    state::AbstractArray,
    from_state::Type{MEE},
    to_state::Type{Cartesian},
    mu::Real,
)
    length(state) == 6 || throw(ArgumentError("State must be of length 6."))
    return _mee_to_cart_partials(state, mu)
end

# Define low-level conversion functions
function _cart_to_kep(state, mu)
    # Position and velocity norms
    r = @inbounds sqrt(state[1]*state[1] + state[2]*state[2] + state[3]*state[3])
    v = @inbounds sqrt(state[4]*state[4] + state[5]*state[5] + state[6]*state[6])

    # Precompute repeated computations
    rinv   = 1.0 / r
    muinv  = 1.0 / mu
    vsq    = v*v
    rDotV  = @inbounds state[1]*state[4] + state[2]*state[5] + state[3]*state[6]
    muRinv = mu*rinv
    twoπ   = 2.0*pi

    # Special case tolerance
    scTol = 1.0e-14;

    # Angular momentum
    hv = @inbounds SVector(
        state[2]*state[6] - state[3]*state[5],
        state[3]*state[4] - state[1]*state[6],
        state[1]*state[5] - state[2]*state[4],
    )
    h  = @inbounds sqrt(hv[1]*hv[1] + hv[2]*hv[2] + hv[3]*hv[3])

    # Normal vector
    nv = @inbounds SVector(-hv[2], hv[1], 0.0)
    n  = @inbounds sqrt(nv[1]^2 + nv[2]^2)

    # Eccentricity
    term1 = (vsq - muRinv)*muinv
    term2 = rDotV*muinv
    ev = @inbounds SVector(
        term1*state[1] - term2*state[4],
        term1*state[2] - term2*state[5],
        term1*state[3] - term2*state[6],
    )
    e = @inbounds sqrt(ev[1]*ev[1] + ev[2]*ev[2] + ev[3]*ev[3])

    # Specific energy
    enrgy = 0.5*vsq - muRinv

    # Semiparameter
    #flag = e != 1.0
    #a = flag ? -mu / (2.0*enrgy) : Inf
    a = -mu / (2.0*enrgy)
    #p = flag ? a*(1.0 - e*e) : h*h*muinv

    # Inclination
    i = @inbounds acos(hv[3] / h)

    # Check if elliptical
    if e > scTol

        # RAAN
        Ω = @inbounds (nv[2] >= 0.0) ? acos(nv[1] / n) : twoπ - acos(nv[1] / n)

        # True anomoly
        eDotR = @inbounds ev[1]*state[1] + ev[2]*state[2] + ev[3]*state[3]
        term  = eDotR / (r*e)
        #if abs(term) > 1.0
        #    term = sign(term)
        #end
        ν = (rDotV >= 0.0) ? acos(term) : twoπ - acos(term)

        # Argument of periapsis/True longitude of periapsis
        if i > scTol && abs(i - π) > scTol
            nDotE = @inbounds nv[1]*ev[1] + nv[2]*ev[2] + nv[3]*ev[3]
            tt    = nDotE / (n*e)
            #tt    = (tt >  1.0 ? 1.0 : tt)
            #tt    = (tt < -1.0 ? -1.0 : tt)
            ω     = @inbounds (ev[3] >= 0.0) ? acos(tt) : (twoπ - acos(tt))

            return SVector(a, e, i, Ω, ω, ν)

        else # Elliptical equatorial
            # For now, will throw an error for singular cases in the Keplerian state
            # representation. We can change this begaviour later if necessary
            throw(DomainError("state corresponds to a equatorial orbit and therefore" *
                "cannot be converted to a Keplerian state."))
        end
    else # Circular
        # For now, will throw error for singular cases in the Keplerian state representation.
        # We can change this behaviour later if necessary
        throw(DomainError("state corresponds to a circular orbit and therefore cannot be" *
            "converted to a Keplerian state."))
    end
end

function _kep_to_cart(state, mu)
    # Compute semi-parameter
    p       = @inbounds state[1]*(1.0 - state[2]^2)

    # Compute position and velocity in PQW frame
    term1   = @inbounds 1.0 / (1.0 + state[2]*cos(state[6]))
    term2   = sqrt(mu / p)
    rpqw    = @inbounds SVector(p*cos(state[6])*term1, p*sin(state[6])*term1, 0.0)
    vpqw    = @inbounds SVector(-term2*sin(state[6]), term2*(state[2] + cos(state[6])), 0.0)

    # Rotate to inertial reference frame
    r       = @inbounds _rot_3_vec(
        _rot_1_vec(_rot_3_vec(rpqw, -state[5]), -state[3]), -state[4]
    )
    v       = @inbounds _rot_3_vec(
        _rot_1_vec(_rot_3_vec(vpqw, -state[5]), -state[3]), -state[4]
    )
    return @inbounds SVector(r[1], r[2], r[3], v[1], v[2], v[3])
end

function _cart_to_mee(state, mu)
    # Compute requirements
    r       = norm(view(state, 1:3))
    v       = norm(view(state, 4:6))
    rdotv   = dot(view(state,1:3),view(state,4:6))
    rhat    = SVector(state[1] / r, state[2] / r, state[3] / r)

    hvec    = @inbounds SVector(
        state[2]*state[6] - state[3]*state[5],
        state[3]*state[4] - state[1]*state[6],
        state[1]*state[5] - state[2]*state[4],
    )
    hmag    = @inbounds sqrt(hvec[1]*hvec[1] + hvec[2]*hvec[2] + hvec[3]*hvec[3])
    hhat    = @inbounds SVector(hvec[1] / hmag, hvec[2] / hmag, hvec[3] / hmag)
    vhat    = @inbounds SVector(
        (r*state[4] - rdotv*rhat[1]) / hmag,
        (r*state[5] - rdotv*rhat[2]) / hmag,
        (r*state[6] - rdotv*rhat[3]) / hmag,
    )
    muInv   = 1.0 / mu

    # Compute elements
    p       = hmag*hmag * muInv
    k       = @inbounds  hhat[1] / (1.0 + hhat[3])
    h       = @inbounds -hhat[2] / (1.0 + hhat[3])

    kk      = k*k
    hh      = h*h
    s2      = 1.0 + hh + kk
    tkh     = 2.0*k*h

    ecc     = @inbounds SVector(
        (state[5]*hvec[3] - state[6]*hvec[2]) * muInv - rhat[1],
        (state[6]*hvec[1] - state[4]*hvec[3]) * muInv - rhat[2],
        (state[4]*hvec[2] - state[5]*hvec[1]) * muInv - rhat[3]
    )

    s2Inv   = 1.0 / s2
    fhat1   = (1.0 - kk + hh) * s2Inv
    fhat2   = tkh * s2Inv
    fhat3   = (-2.0 * k) * s2Inv
    ghat1   = tkh * s2Inv
    ghat2   = (1.0 + kk - hh) * s2Inv
    ghat3   = (2.0 * h) * s2Inv

    f       = @inbounds ecc[1]*fhat1 + ecc[2]*fhat2 + ecc[3]*fhat3
    g       = @inbounds ecc[1]*ghat1 + ecc[2]*ghat2 + ecc[3]*ghat3

    sinl    = rhat[2] - vhat[1]
    cosl    = rhat[1] + vhat[2]
    L       = atan(sinl, cosl)

    return SVector(p,f,g,h,k,L)
end

function _mee_to_cart(state, mu)
    α2      = @inbounds state[4]*state[4] - state[5]*state[5]
    s2      = @inbounds 1.0 + state[4]*state[4] + state[5]*state[5]
    s2Inv   = 1.0 / s2
    w       = @inbounds 1.0 + state[2]*cos(state[6]) + state[3]*sin(state[6])
    r       = @inbounds state[1] / w
    sqrtmup = @inbounds sqrt(mu / state[1])

    cart    = @inbounds SVector(
        r*s2Inv*(cos(state[6]) + α2*cos(state[6]) + 2*state[4]*state[5]*sin(state[6])),
        r*s2Inv*(sin(state[6]) - α2*sin(state[6]) + 2*state[4]*state[5]*cos(state[6])),
        2.0*r*s2Inv*(state[4]*sin(state[6]) - state[5]*cos(state[6])),
        -s2Inv*sqrtmup*(
            sin(state[6]) + α2*sin(state[6]) - 2*state[4]*state[5]*cos(state[6]) +
            state[3] - 2.0*state[2]*state[4]*state[5] + α2*state[3]
        ),
        -s2Inv*sqrtmup*(
            -cos(state[6]) + α2*cos(state[6]) + 2*state[4]*state[5]*sin(state[6]) -
            state[2] + 2.0*state[3]*state[4]*state[5] + α2*state[2]
        ),
        2.0*s2Inv*sqrtmup*(
            state[4]*cos(state[6]) + state[5]*sin(state[6]) +
            state[2]*state[4] + state[3]*state[5]
        ),
    )
    return cart
end

function _kep_to_mee(state, mu)
    p = @inbounds state[1]*(1.0 - state[2]*state[2])
    f = @inbounds state[2]*cos(state[4] + state[5])
    g = @inbounds state[2]*sin(state[4] + state[5])
    h = @inbounds tan(0.5*state[3])*cos(state[4])
    k = @inbounds tan(0.5*state[3])*sin(state[4])
    L = @inbounds state[4] + state[5] + state[6]
    return SVector(p,f,g,h,k,L)
end

function _mee_to_kep(state, mu)
    a = @inbounds state[1] / (1.0 - state[2]*state[2] - state[3]*state[3])
    e = @inbounds sqrt(state[2]*state[2] + state[3]*state[3])
    i = @inbounds atan(
        2.0*sqrt(state[4]*state[4] + state[5]*state[5]),
        1.0 - state[4]*state[4] - state[5]*state[5],
    )
    Ω = @inbounds atan(state[5], state[4])
    ω = @inbounds atan(
        state[3]*state[4] - state[2]*state[5],
        state[2]*state[4] + state[3]*state[5],
    )
    ν = @inbounds state[6] - atan(state[3] / state[2])
    return SVector(a,e,i,Ω,ω,ν)
end

function _mee_to_cart_partials(state, mu)
    # Instantiate variables for generated code
    p,f,g,h,k,L = state

    # Start of generated code
    t2 = cos(L);
    t3 = sin(L);
    t4 = f*h;
    t5 = g*k;
    t6 = h*h;
    t7 = h*h*h;
    t8 = k*k;
    t9 = k*k*k;
    t18 = 1.0/p;
    t10 = f*t2;
    t11 = g*t2;
    t12 = h*t2;
    t13 = k*t2;
    t14 = f*t3;
    t15 = g*t3;
    t16 = h*t3;
    t17 = k*t3;
    t19 = t18*t18;
    t20 = -t2;
    t21 = k*t4*2.0;
    t22 = h*t5*2.0;
    t23 = -t8;
    t27 = mu*t18;
    t28 = t2*t6;
    t29 = t2*t8;
    t31 = t3*t6;
    t32 = t3*t8;
    t37 = t6+t8+1.0;
    t24 = -t14;
    t25 = -t16;
    t26 = -t21;
    t30 = k*t12*2.0;
    t33 = k*t16*2.0;
    t35 = t8*t20;
    t36 = -t31;
    t38 = t6+t23;
    t39 = t10+t15+1.0;
    t40 = sqrt(t27);
    t41 = 1.0/t37;
    t51 = t4+t5+t12+t17;
    t34 = -t30;
    t42 = t41*t41;
    t43 = t11+t24;
    t44 = t13+t25;
    t45 = 1.0/t40;
    t46 = t2*t38;
    t47 = t3*t38;
    t48 = 1.0/t39;
    t53 = t2+t28+t33+t35;
    t54 = t3+t30+t32+t36;
    t56 = h*k*t40*t41*2.0;
    t49 = t48*t48;
    t50 = -t47;
    t52 = t2+t33+t46;
    t55 = t3+t30+t50;

    A11 = t41*t48*t52;
    A12 = p*t20*t41*t49*t52;
    A13 = -p*t3*t41*t49*t52;
    A14 = k*p*t42*t48*t54*2.0;
    A15 = p*t42*t48*(t13+t44-t3*t7+t6*t13*2.0+t8*t16)*-2.0;
    A16 = -p*t41*t48*(t3+t34+t47)-p*t41*t43*t49*t52;
    A21 = t41*t48*t55;
    A22 = p*t20*t41*t49*t55;
    A23 = -p*t3*t41*t49*t55;
    A24 = p*t42*t48*(-t13+t16*2.0+t6*t13+t8*t16*2.0+t9*t20)*-2.0;
    A25 = h*p*t42*t48*t53*2.0;
    A26 = -p*t41*t48*(t20+t33+t46)-p*t41*t43*t49*t55;
    A31 = t41*t44*t48*-2.0;
    A32 = p*t2*t41*t44*t49*2.0;
    A33 = p*t3*t41*t44*t49*2.0;
    A34 = p*t42*t48*t54*2.0;
    A35 = p*t42*t48*t53*-2.0;
    A36 = p*t41*t49*t51*2.0;
    A41 = (mu*t19*t41*t45*(g+t3+t26+t34+t47+g*t38))/2.0;
    A42 = t56;
    A43 = -t40*t41*(t38+1.0);
    A44 = k*t40*t42*(f+t2-t22+t29-t33+f*t8-h*t4+t6*t20)*2.0;
    A45 = t40*t42*(t5+t17+t51+t2*t7+t4*t6+t5*t6*2.0+t6*t17*2.0+t4*t23+t12*t23)*2.0;
    A46 = -t40*t41*t52;
    A51 = (mu*t19*t41*t45*(-f+t20+t22+t33+t46+f*t38))/2.0;
    A52 = t40*t41*(-t6+t8+1.0);
    A53 = -t56;
    A54 = t40*t42*(t4+t12+t51-t5*t6+t3*t9+t4*t8*2.0+t5*t8+t8*t12*2.0-t6*t17)*-2.0;
    A55 = h*t40*t42*(g+t3+t26+t31+t34+g*t6-k*t5+t3*t23)*-2.0;
    A56 = -t40*t41*t55;
    A61 = -mu*t19*t41*t45*t51;
    A62 = h*t40*t41*2.0;
    A63 = k*t40*t41*2.0;
    A64 = t40*t41*(f+t2)*2.0-h*t40*t42*t51*4.0;
    A65 = t40*t41*(g+t3)*2.0-k*t40*t42*t51*4.0;
    A66 = t40*t41*t44*2.0;
    return SMatrix{6,6}(A11,A21,A31,A41,A51,A61,
                        A12,A22,A32,A42,A52,A62,
                        A13,A23,A33,A43,A53,A63,
                        A14,A24,A34,A44,A54,A64,
                        A15,A25,A35,A45,A55,A65,
                        A16,A26,A36,A46,A56,A66)
end
