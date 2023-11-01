
# ===== Linear solves
"""
    _linsolve(
        x::AbstractVector,
        A::AbstractMatrix,
        b::AbstractVector,
    )

Solves the linear system Ax = b for x using a pure Julia implementation of an LU factorization.

# Arguments
- `x::AbstractVector`: Vector to store solution in.
- `A::AbstractMatrix`: The A matrix of the linear system.
- `b::AbstractVector`: The b vector of the linear system.
"""
function _linsolve!(x, A, b)
    # Preconditioners
    Pl = I
    Pr = I

    # Solve the linear system
    p   = LinearProblem(A, b)
    s   = solve(p, RFLUFactorization(); Pl = Pl, Pr = Pr)
    x  .= s.u
    return nothing
end

"""
    _linsolve(
        x::AbstractVector,
        A::SparseMatrixCSC,
        b::AbstractVector,
    )

Solves the linear system Ax = b for x using a Krylov GMRES iterative solver.

# Arguments
- `x::AbstractVector`: Vector to store solution in.
- `A::AbstractMatrix`: The A matrix of the linear system.
- `b::AbstractVector`: The b vector of the linear system.
"""
function _linsolve!(x, A::SparseMatrixCSC, b)
    # Preconditioners
    Pl = I
    Pr = I

    # Solve the linear system
    p   = LinearProblem(A, b)
    s   = solve(p, KrylovJL_GMRES(); Pl = Pl, Pr = Pr)
    x  .= s.u
    return nothing
end

# ===== Vector rotations

"""
    _rot_1_vec(
        vec_in::AbstractVector,
        a::Real,
    )

Rotates `vec_in` by angle `a` about 1st axis.

# Arguments
- `vec_in::AbstractVector`: Vector to rotate.
- `a::Real`: Angle to rotate by (radians).

# Returns
- `vec_out::AbstractVector`: Rotated vector.
"""
function _rot_1_vec(vec_in, a)
    c = cos(a)
    s = sin(a)
    @inbounds vec_out = SVector(
        vec_in[1],
        c*vec_in[2] + s*vec_in[3],
        c*vec_in[3] - s*vec_in[2],
    )
    return vec_out
end

"""
    _rot_1_vec!(
        vec::AbstractVector,
        a::Real,
    )

Rotates `vec` in-place by angle `a` about 1st axis.

# Arguments
- `vec::AbstractVector`: Vector to rotate.
- `a::Real`: Angle to rotate by (radians).

# Returns
- `nothing`
"""
function _rot_1_vec!(vec, a)
    c  = cos(a)
    s  = sin(a)
    v2 = @inbounds c*vec_in[2] + s*vec_in[3]
    v3 = @inbounds c*vec_in[3] - s*vec_in[2]
    @inbounds vec[2] = v2
    @inbounds vec[3] = v3
    return nothing
end

"""
    _rot_2_vec(
        vec_in::AbstractVector,
        a::Real,
    )

Rotates `vec_in` by angle `a` about 2nd axis.

# Arguments
- `vec_in::AbstractVector`: Vector to rotate.
- `a::Real`: Angle to rotate by (radians).

# Returns
- `vec_out::AbstractVector`: Rotated vector.
"""
function _rot_2_vec(vec_in, a)
    c = cos(a)
    s = sin(a)
    @inbounds vec_out = SVector(
        c*vec_in[1] - s*vec_in[3],
        vec_in[2],
        c*vec_in[3] + s*vec_in[1],
    )
    return vec_out
end

"""
    _rot_2_vec!(
        vec::AbstractVector,
        a::Real,
    )

Rotates `vec` in-place by angle `a` about 2nd axis.

# Arguments
- `vec::AbstractVector`: Vector to rotate.
- `a::Real`: Angle to rotate by (radians).

# Returns
- `nothing`
"""
function _rot_2_vec!(vec, a)
    c  = cos(a)
    s  = sin(a)
    v1 = c*vec_in[1] - s*vec_in[3]
    v3 = c*vec_in[3] + s*vec_in[1]
    vec[1] = v1
    vec[3] = v3
    return nothing
end

"""
    _rot_3_vec(
        vec_in::AbstractVector,
        a::Real,
    )

Rotates `vec_in` by angle `a` about 3rd axis.

# Arguments
- `vec_in::AbstractVector`: Vector to rotate.
- `a::Real`: Angle to rotate by (radians).

# Returns
- `vec_out::AbstractVector`: Rotated vector.
"""
function _rot_3_vec(vec_in, a)
    c = cos(a)
    s = sin(a)
    @inbounds vec_out = SVector(
        c*vec_in[1] + s*vec_in[2],
        c*vec_in[2] - s*vec_in[1],
        vec_in[3],
    )
    return vec_out
end

"""
    _rot_3_vec!(
        vec::AbstractVector,
        a::Real,
    )

Rotates `vec` in-place by angle `a` about 3rd axis.

# Arguments
- `vec::AbstractVector`: Vector to rotate.
- `a::Real`: Angle to rotate by (radians).

# Returns
- `nothing``
"""
function _rot_3_vec!(vec, a)
    c  = cos(a)
    s  = sin(a)
    v1 = c*vec_in[1] + s*vec_in[2]
    v2 = c*vec_in[2] - s*vec_in[1]
    vec[1] = v1
    vec[2] = v2
    return nothing
end

"""
    _inertial_to_rtn_vec(
        vec::AbstractVector,
        cart_state::AbstractVector,
    )

Rotates `vec` from an inertial frame to the local RTN frame of a spacecraft, given its
cartesian state vector `cart_state`, defined w.r.t. the inertial frame.

# Arguments
- `vec::AbstractVector`: Vector to rotate.
- `cart_state::AbstractVector`: Cartesian state vector of spacecraft [pos; vel].

# Returns
- `vec_out::AbstractVector`: Rotated vector.
"""
function _inertial_to_rtn_vec(vec, cart_state)
    # Get position and velocity from cartesian state
    rvec = view(cart_state, 1:3)
    vvec = view(cart_state, 4:6)

    # Construct rhat vector
    rinv = 1.0 / norm(rvec)
    rhat = @inbounds SVector(rinv*rvec[1], rinv*rvec[2], rinv*rvec[3])

    # Construct hhat vector
    hvec = cross(rvec, vvec)
    hinv = 1.0 / norm(hvec)
    nhat = @inbounds SVector(hinv*hvec[1], hinv*hvec[2], hinv*hvec[3])

    # Construct that vector
    tvec = cross(nhat, rhat)
    tinv = 1.0 / norm(tvec)
    that = @inbounds SVector(tinv*tvec[1], tinv*tvec[2], tinv*tvec[3])

    # Compute accelerations by projecting onto each basis vector of RTN frame
    return SVector(
        dot(rhat, vec),
        dot(that, vec),
        dot(nhat, vec),
    )
end

"""
    _inertial_to_rtn_vec_partial_cart(
        vec::AbstractVector,
        dvecdcart::AbstractMatrix,
        cart_state::AbstractVector,
    )

Computes the partial derivative of R(`cart_state`) * `vec`(`cart_state`) w.r.t
`cart_state`, where R is the rotation matrix from an inertial frame to the local RTN
frame of a spacecraft.

# Arguments
- `vec::AbstractVector`: Vector to rotate.
- `dvecdcart::AbstractMatrix`: Partial derivative of `vec` w.r.t. `cart_state`.
- `cart_state::AbstractVector`: Cartesian state vector of spacecraft [pos; vel].

# Returns
- `out::StaticMatrix{3,6,Float64}`: The partial derivative.
"""
function _inertial_to_rtn_vec_partial_cart(vec, dvecdcart, cart_state)
    # Get requirements for generated code
    x1,x2,x3,x4,x5,x6 = cart_state
    vec1,vec2,vec3 = vec
    dvecdx1_1,dvecdx1_2,dvecdx1_3,dvecdx1_4,dvecdx1_5,dvecdx1_6 = view(dvecdcart, 1,:)
    dvecdx2_1,dvecdx2_2,dvecdx2_3,dvecdx2_4,dvecdx2_5,dvecdx2_6 = view(dvecdcart, 2,:)
    dvecdx3_1,dvecdx3_2,dvecdx3_3,dvecdx3_4,dvecdx3_5,dvecdx3_6 = view(dvecdcart, 3,:)

    # Start of generated code
    t2 = x1*x5;
    t3 = x2*x4;
    t4 = x1*x6;
    t5 = x3*x4;
    t6 = x2*x6;
    t7 = x3*x5;
    t8 = x1*x1;
    t9 = x1*x1*x1;
    t10 = x2*x2;
    t12 = x2*x2*x2;
    t13 = x3*x3;
    t16 = x3*x3*x3;
    t24 = x6*x6;
    t19 = t13*t13;
    t26 = t6*vec1;
    t27 = t7*vec1;
    t28 = t4*vec2;
    t29 = t5*vec2;
    t30 = t2*vec3;
    t31 = t3*vec3;
    t32 = t3*x1;
    t33 = t2*x2;
    t34 = t5*x1;
    t35 = t4*x3;
    t36 = t7*x2;
    t37 = t6*x3;
    t42 = t3*vec1*x5;
    t43 = t2*vec2*x4;
    t46 = t5*vec1*x6;
    t47 = t4*vec3*x4;
    t48 = t7*vec2*x6;
    t49 = t6*vec3*x5;
    t50 = -t3;
    t51 = -t5;
    t52 = -t7;
    t53 = t2*x1;
    t54 = t3*x2;
    t55 = t4*x1;
    t56 = t5*x3;
    t57 = t6*x2;
    t58 = t7*x3;
    t59 = t2*vec1*x5;
    t61 = t4*vec1*x6;
    t65 = t6*vec2*x6;
    t71 = t2*t3*2.0;
    t72 = t4*t5*2.0;
    t73 = t6*t7*2.0;
    t83 = t2*t2;
    t84 = t3*t3;
    t85 = t4*t4;
    t86 = t5*t5;
    t87 = t6*t6;
    t88 = t7*t7;
    t107 = t8+t10+t13;
    t38 = t33*vec1;
    t39 = t32*vec2;
    t40 = t35*vec1;
    t41 = t34*vec3;
    t44 = t37*vec2;
    t45 = t36*vec3;
    t60 = t54*vec1;
    t62 = t56*vec1;
    t64 = t53*vec2;
    t66 = t58*vec2;
    t68 = t55*vec3;
    t70 = t57*vec3;
    t74 = -t27;
    t75 = -t28;
    t76 = -t31;
    t77 = -t32;
    t78 = -t33;
    t82 = -t37;
    t101 = -t71;
    t102 = -t72;
    t103 = -t73;
    t104 = t2+t50;
    t105 = t4+t51;
    t106 = t6+t52;
    t111 = sqrt(t107);
    t95 = -t44;
    t108 = t104*t104;
    t109 = t105*t105;
    t110 = t106*t106;
    t112 = 1.0/t111;
    t117 = t26+t29+t30+t74+t75+t76;
    t120 = t83+t84+t85+t86+t87+t88+t101+t102+t103;
    t113 = t112*t112*t112;
    t118 = t108+t109+t110;
    t121 = 1.0/sqrt(t120);
    t119 = 1.0/sqrt(t118);
    t122 = t121*t121*t121;
    A011 = t113*(t10*vec1+t13*vec1+dvecdx1_1*t107*x1+dvecdx2_1*t107*x2+dvecdx3_1*t107*x3-vec2*x1*x2-vec3*x1*x3);
    A012 = t113*(t8*vec2+t13*vec2+dvecdx1_2*t107*x1+dvecdx2_2*t107*x2+dvecdx3_2*t107*x3-vec1*x1*x2-vec3*x2*x3);
    A013 = t113*(t8*vec3+t10*vec3+dvecdx1_3*t107*x1+dvecdx2_3*t107*x2+dvecdx3_3*t107*x3-vec1*x1*x3-vec2*x2*x3);
    A014 = t112*(dvecdx1_4*x1+dvecdx2_4*x2+dvecdx3_4*x3);
    A015 = t112*(dvecdx1_5*x1+dvecdx2_5*x2+dvecdx3_5*x3);
    A016 = t112*(dvecdx1_6*x1+dvecdx2_6*x2+dvecdx3_6*x3);
    A021 = t113*t122*((t2*t2*t2)*t10*vec2-(t3*t3*t3)*t13*vec2*2.0+(t4*t4*t4)*t13*vec3-t2*t32*t38*3.0+t3*t32*t38*3.0-t2*t32*t40*2.0+t3*t32*t40*3.0-t4*t34*t40*3.0+t5*t34*t40*3.0+t2*t38*t53+t2*t40*t53+t2*t32*t62*3.0-t3*t32*t62*2.0+t4*t38*t55+t4*t40*t55+t26*t36*t37-t30*t36*t37*2.0+t31*t36*t37*2.0+t26*t36*t58+t28*t37*t57*2.0+t30*t36*t58*2.0-t31*t36*t58*2.0+t32*t50*t60+t34*t51*t62+t26*t57*t82+t3*t4*t8*t28-t3*t4*t10*t26*2.0+t2*t4*t10*t28*2.0+t2*t5*t8*t29+t2*t5*t8*t30-t2*t4*t10*t30+t3*t4*t8*t31-t3*t4*t13*t26*2.0-t2*t5*t13*t27*2.0+t3*t4*t10*t31-t3*t4*t13*t28*3.0+t2*t4*t13*t30*2.0+t2*t5*t13*t29+t3*t5*t13*t28*2.0-t3*t4*t13*t30*6.0+t3*t4*t13*t31*4.0+t3*t5*t13*t30*2.0-t3*t5*t13*t31*2.0+t2*t8*t28*t51+t5*t8*t28*t50+t4*t8*t30*t50+t5*t8*t30*t50+t2*t4*t13*t75+t5*t13*t29*t50-t2*t19*t24*vec2+t3*t13*t73*vec2+t19*t24*t50*vec2+t2*t10*t84*vec2*3.0-t3*t10*t83*vec2*3.0+t2*t10*t87*vec2*2.0+t2*t13*t84*vec2*4.0+t2*t10*t88*vec2-t2*t13*t87*vec2*2.0-t3*t13*t87*vec2*2.0+t4*t13*t86*vec3*3.0-t5*t13*t85*vec3*3.0+t4*t13*t87*vec3+t10*t50*t84*vec2+t10*t50*t87*vec2+t13*t51*t86*vec3+t13*t51*t88*vec3-t6*t12*t26*x5+t7*t12*t26*x5*2.0-t6*t12*t30*x5+t6*t12*t31*x5-t6*t16*t26*x6+t7*t12*t30*x5+t7*t16*t26*x6*2.0+t6*t16*t28*x6-t6*t16*t30*x6*3.0+t7*t16*t29*x6+t7*t16*t30*x6*2.0+t7*t16*t31*x6+t3*t26*t33*x3*4.0-t3*t27*t33*x3*2.0-t3*t28*t33*x3*6.0+t2*t30*t33*x3*2.0-t3*t30*t33*x3*3.0+t4*t28*t35*x2*2.0+t3*t28*t54*x3*2.0+t3*t30*t54*x3*2.0+t12*t27*t52*x5+t12*t31*t52*x5+t16*t27*t52*x6+t27*t36*t52*x3+t31*t50*t54*x3-t32*t85*vec1*x2*2.0-t34*t83*vec1*x3*2.0-t2*t6*t7*t10*vec2*3.0+t3*t6*t7*t10*vec2+t2*t3*t16*vec1*x6*4.0)-dvecdx2_1*t112*t121*(t32+t37-t53+t52*x3)-dvecdx1_1*t112*t121*(t33+t35+t50*x2+t51*x3)-dvecdx3_1*t112*t121*(t34+t36-t55-t57);
    A022 = -t113*t122*((t2*t2*t2)*t8*vec1+(t2*t2*t2)*t13*vec1*2.0-(t6*t6*t6)*t13*vec3+(t7*t7*t7)*t13*vec3-t2*t33*t39*3.0-t2*t35*t38*2.0+t3*t33*t39*3.0+t3*t35*t38*6.0-t2*t33*t44*3.0+t3*t33*t44*2.0+t6*t36*t44*3.0-t7*t36*t44*3.0+t2*t33*t64-t30*t34*t35*2.0+t31*t34*t35*2.0+t2*t33*t66*2.0-t3*t33*t66*3.0+t7*t36*t66+t28*t33*t55*2.0-t26*t35*t57*2.0+t28*t35*t55+t29*t34*t56+t30*t34*t56*2.0-t31*t34*t56*2.0+t39*t50*t54+t34*t35*t75+t44*t50*t54+t6*t57*t95+t30*t55*t78-t40*t85*x2*2.0-t2*t6*t10*t26+t2*t7*t10*t26+t2*t6*t10*t28*2.0+t3*t7*t10*t26+t2*t6*t13*t26*3.0-t2*t6*t10*t30-t2*t7*t13*t26*2.0+t3*t6*t13*t26+t2*t6*t13*t28*2.0+t2*t7*t13*t27+t3*t6*t10*t30+t3*t7*t10*t30-t2*t6*t13*t30*4.0+t2*t7*t13*t30*2.0+t3*t6*t13*t30*6.0+t3*t7*t13*t29*2.0-t3*t6*t13*t31*2.0-t3*t7*t13*t30*2.0+t6*t10*t28*t50+t7*t10*t27*t50+t7*t13*t27*t50+t7*t10*t31*t50+t2*t19*t24*vec1+t3*t19*t24*vec1+t2*t8*t84*vec1*3.0-t3*t8*t83*vec1*3.0+t2*t8*t85*vec1-t3*t8*t85*vec1*2.0-t3*t10*t85*vec1*2.0-t3*t13*t83*vec1*4.0+t2*t13*t85*vec1*2.0+t3*t13*t85*vec1*2.0-t6*t13*t85*vec3+t7*t13*t86*vec3-t6*t13*t88*vec3*3.0+t7*t13*t87*vec3*3.0+t8*t50*t84*vec1+t8*t50*t86*vec1+t4*t9*t28*x4-t5*t9*t28*x4*2.0-t4*t9*t30*x4+t5*t9*t29*x4+t4*t9*t31*x4+t5*t9*t30*x4-t4*t16*t26*x6+t4*t16*t28*x6-t5*t16*t28*x6*2.0+t5*t16*t29*x6+t4*t16*t31*x6*3.0-t5*t16*t31*x6*2.0-t2*t28*t32*x3*4.0+t2*t29*t32*x3*2.0-t2*t30*t32*x3*2.0+t3*t30*t32*x3*3.0-t3*t31*t32*x3*2.0+t4*t31*t32*x2+t2*t30*t53*x3+t9*t31*t51*x4+t16*t27*t51*x6+t16*t30*t51*x6+t28*t34*t51*x3+t36*t84*vec2*x3*2.0+t3*t4*t5*t8*vec1*3.0-t2*t4*t5*t13*vec1*2.0+t2*t4*t8*t51*vec1-t2*t3*t16*vec2*x6*4.0)-dvecdx2_2*t112*t121*(t32+t37-t53+t52*x3)-dvecdx1_2*t112*t121*(t33+t35+t50*x2+t51*x3)-dvecdx3_2*t112*t121*(t34+t36-t55-t57);
    A023 = -t113*t122*((t4*t4*t4)*t8*vec1+(t4*t4*t4)*t10*vec1*2.0+(t6*t6*t6)*t10*vec2-t4*t35*t38*2.0-t4*t35*t41*3.0+t5*t35*t41*3.0-t6*t37*t45*3.0+t7*t37*t45*3.0-t28*t32*t33*2.0-t28*t32*t35*2.0-t26*t33*t37*2.0-t30*t32*t35*4.0-t4*t40*t54*4.0-t28*t33*t37*4.0+t31*t32*t35*2.0+t30*t33*t37*2.0+t4*t35*t68+t4*t35*t70*2.0+t6*t37*t70+t28*t32*t54*2.0+t28*t32*t56*3.0+t30*t33*t53+t26*t33*t58*3.0-t29*t32*t56*2.0+t31*t32*t54-t27*t33*t58*2.0+t30*t35*t53*2.0-t28*t37*t54*2.0-t30*t37*t54*4.0+t31*t37*t54*2.0+t30*t33*t77+t41*t51*t56+t45*t52*t58+t35*t53*t75+t38*t71*x3+t39*t71*x3-t38*t83*x3*2.0-t39*t84*x3*2.0+t4*t6*t10*t26+t4*t6*t10*t28*2.0-t4*t6*t13*t30*3.0-t4*t6*t13*t31*3.0+t7*t13*t30*t51+t7*t13*t31*t51+t4*t8*t83*vec1-t5*t8*t83*vec1*2.0+t4*t10*t83*vec1*2.0+t4*t8*t86*vec1*3.0-t5*t8*t85*vec1*3.0-t5*t13*t83*vec1*2.0+t6*t10*t84*vec2-t7*t10*t84*vec2*2.0+t6*t10*t88*vec2*3.0-t7*t10*t87*vec2*3.0-t7*t13*t84*vec2*2.0+t8*t51*t86*vec1+t10*t52*t83*vec2+t10*t52*t88*vec2+t2*t9*t29*x4+t3*t9*t28*x4+t2*t9*t30*x4+t2*t12*t26*x5-t3*t9*t30*x4*2.0+t3*t9*t31*x4+t3*t12*t27*x5+t2*t12*t30*x5+t2*t16*t26*x6-t3*t12*t30*x5*2.0+t3*t12*t31*x5+t3*t16*t27*x6+t2*t16*t29*x6+t3*t16*t28*x6+t2*t16*t30*x6*2.0+t3*t16*t30*x6*2.0+t3*t16*t31*x6*2.0+t2*t29*t34*x3+t3*t27*t36*x3+t4*t28*t55*x2+t2*t9*t75*x4+t9*t29*t50*x4+t2*t12*t74*x5+t12*t26*t50*x5+t2*t16*t74*x6+t16*t26*t50*x6+t2*t16*t75*x6+t16*t29*t50*x6+t30*t32*t50*x2+t26*t37*t50*x2+t2*t3*t5*t8*vec1*3.0-t2*t3*t4*t10*vec1*2.0+t2*t3*t4*t13*vec1*6.0+t2*t3*t7*t10*vec2*3.0+t2*t3*t6*t13*vec2*6.0+t2*t4*t8*t50*vec1+t3*t5*t8*t50*vec1+t2*t6*t10*t50*vec2)-dvecdx2_3*t112*t121*(t32+t37-t53+t52*x3)-dvecdx1_3*t112*t121*(t33+t35+t50*x2+t51*x3)-dvecdx3_3*t112*t121*(t34+t36-t55-t57);
    A024 = -dvecdx2_4*t112*t121*(t32+t37-t53+t52*x3)+t106*t111*t117*t122-dvecdx1_4*t112*t121*(t33+t35+t50*x2+t51*x3)-dvecdx3_4*t112*t121*(t34+t36-t55-t57);
    A025 = -dvecdx2_5*t112*t121*(t32+t37-t53+t52*x3)-t105*t111*t117*t122-dvecdx1_5*t112*t121*(t33+t35+t50*x2+t51*x3)-dvecdx3_5*t112*t121*(t34+t36-t55-t57);
    A026 = -dvecdx2_6*t112*t121*(t32+t37-t53+t52*x3)+t104*t111*t117*t122-dvecdx1_6*t112*t121*(t33+t35+t50*x2+t51*x3)-dvecdx3_6*t112*t121*(t34+t36-t55-t57);
    A031 = t106*t122*(t42+t43+t46+t47+t48+t49-t59-t61-t65+t50*vec2*x4+t51*vec3*x4+t52*vec3*x5)+dvecdx1_1*t106*t119-dvecdx2_1*t105*t119+dvecdx3_1*t104*t119;
    A032 = -t105*t122*(t42+t43+t46+t47+t48+t49-t59-t61-t65+t50*vec2*x4+t51*vec3*x4+t52*vec3*x5)+dvecdx1_2*t106*t119-dvecdx2_2*t105*t119+dvecdx3_2*t104*t119;
    A033 = t104*t122*(t42+t43+t46+t47+t48+t49-t59-t61-t65+t50*vec2*x4+t51*vec3*x4+t52*vec3*x5)+dvecdx1_3*t106*t119-dvecdx2_3*t105*t119+dvecdx3_3*t104*t119;
    A034 = t106*t122*(t38+t39+t40+t41+t44+t45-t64-t68-t70+t50*vec1*x2+t51*vec1*x3+t52*vec2*x3)+dvecdx1_4*t106*t119-dvecdx2_4*t105*t119+dvecdx3_4*t104*t119;
    A035 = -t105*t122*(t38+t39+t40+t41+t44+t45-t64-t68-t70+t50*vec1*x2+t51*vec1*x3+t52*vec2*x3)+dvecdx1_5*t106*t119-dvecdx2_5*t105*t119+dvecdx3_5*t104*t119;
    A036 = t104*t122*(t38+t39+t40+t41+t44+t45-t64-t68-t70+t50*vec1*x2+t51*vec1*x3+t52*vec2*x3)+dvecdx1_6*t106*t119-dvecdx2_6*t105*t119+dvecdx3_6*t104*t119;

    return SMatrix{3,6}(
        A011, A021, A031,
        A012, A022, A032,
        A013, A023, A033,
        A014, A024, A034,
        A015, A025, A035,
        A016, A026, A036,
    )
end
