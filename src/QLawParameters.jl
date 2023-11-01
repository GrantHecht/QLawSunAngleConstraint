# Define some qlaw parameters
mutable struct QLawParameters{DES}
    # Initial orbital elements
    oe0::Vector{Float64}

    # Target orbital elements
    oet::Vector{Float64}

    # Target element weights
    oeW::Vector{Float64}

    # Target element tolerances
    oeTols::Vector{Float64}

    # Target element error convergence switches
    Ws::Vector{Int}

    # QLaw Parameters
    Wp::Float64
    rpmin::Float64
    k::Float64
    ηr::Float64
    ηa::Float64
    steps::Int
    b_petro::Float64
    m_petro::Float64
    n_petro::Float64
    r_petro::Float64

    # Dynamics parameters
    μ::Float64

    # Spacecraft parameters
    m0::Float64
    mp::Float64
    tMax::Float64
    c::Float64

    # QLaw control variables
    α::Float64
    β::Float64
    T::Float64
    coasting::Bool

    # Parameters for the dynamics
    meePs::MEEParams
    scPs::SimpleSpacecraft

    # Thrust-to-sun angle constraint
    thrustSunAngleConstraint::Bool
    toSunVec::Vector{Float64}
    thrustSunAngle::Float64
    panelType::Symbol

    # Solar eclipsing
    eclipsing::Bool
    eclipsed::Bool
    RB::Float64 # Radius of central body
    RS::Float64 # Sun Radius

    # Integration parameters
    desolver::DES
    reltol::Float64
    abstol::Float64
    maxRevs::Float64
    integStep::Float64

    # QLaw type
    type::Symbol

    # Return trajectory at steps
    returnTrajAtSteps::Bool

    # Storage info (Switch to type flag in future)
    writeDataToFile::Bool
    writeDataOnlyAtSteps::Bool
end

function QLawParameters(oe0,oet;
    oeW         = [1.0, 1.0, 1.0, 0.0, 0.0],
    oeTols      = [10.0, 0.001, 0.01, 0.01, 0.01],
    Ws          = nothing,
    rpmin       = 6578.0,
    k           = 100.0,
    Wp          = 0.0,
    b_petro     = 0.01,
    m_petro     = 3.0,
    n_petro     = 4.0,
    r_petro     = 2.0,
    ηr_tol      = 0.0,
    ηa_tol      = 0.0,
    eSteps      = 60,
    meeParams   = MEEParams(0.0, LU = 1.0, MU = 1.0, TU = 1.0),
    spaceCraft  = SimpleSpacecraft(3000.0, 2500.0, 1.0, 3000.0),
    desolver    = Vern9(),
    reltol      = 1e-10,
    abstol      = 1e-10,
    maxRevs     = 100.0,
    integStep   = 5.0,
    returnData  = false,
    writeData   = false,
    onlyWriteDataAtSteps = false,
    type        = :QDUC,
    eclipsing   = false,
    RB          = 6378.14,
    RS          = 695500.0,
    thrustSunAngleConstraint = false,
    thrustSunAngle           = 30.0*pi/180.0,
    toSunVec                 = [1.0, 0.0, 0.0],
    panelType   = :dual # options: :dual, thrustdir, antithrustdir
    )

    # Check argument size
    if length(oe0) != 6
        throw(ArgumentError("Initial orbital element vector is incorrect length."))
    end
    if length(oet) != 5
        throw(ArgumentError("Target orbital element vector is incorrect length."))
    end
    if length(oeW) != 5
        throw(ArgumentError("Target element weights vector is incorrect length."))
    end
    if length(oeTols) != 5 
        throw(ArgumentError("Target orbital element tolerance vector is incorrect length."))
    end
    if !isnothing(Ws) && length(Ws) != 5
        throw(ArgumentError("Target element switch vector is incorrect length."))
    end

    # ===== Scale arguments
    LU  = meeParams.LU
    TU  = meeParams.TU
    MU  = meeParams.MU
    d2r = pi / 180.0

    # Initial orbital elements
    oe0[1]      /= LU
    oe0[3:6]   .*= d2r

    # Target orbital elements
    oet[1]      /= LU
    oet[3:5]   .*= d2r

    # Target element tolerances
    oeTols[1]   /= LU
    oeTols[3:5].*= d2r

    # Other parameters
    rpmin       /= LU
    μ            = AstroEOMs.getScaledGravityParameter(meeParams)
    integStep   *= d2r

    # Spacecraft parameters
    m0           = spaceCraft.initMass / MU
    mp           = spaceCraft.initProp / MU
    tMax         = spaceCraft.tMax * TU*TU / (1000.0*MU*LU)
    c            = 9.80665*spaceCraft.isp * TU / (1000.0 * LU) 

    # Eclipsing parameters
    RB          /= LU
    RS          /= LU

    # Set targeting switches if necessary
    if isnothing(Ws)
        Wa          = oeW[1] > 0.0 ? 1 : 0
        We          = oeW[2] > 0.0 ? 1 : 0
        Wi          = oeW[3] > 0.0 ? 1 : 0
        WΩ          = oeW[4] > 0.0 ? 1 : 0
        Wω          = oeW[5] > 0.0 ? 1 : 0
        Ws          = [Wa,We,Wi,WΩ,Wω]
    end

    # Ensure toSunVec is a unit vector
    mag = norm(toSunVec)
    toSunVec ./= mag

    # Construct parameter type
    QLawParameters(oe0,oet,oeW,oeTols,Ws,Wp,rpmin,k,ηr_tol,ηa_tol,eSteps,
        b_petro,m_petro,n_petro,r_petro,μ,m0,mp,tMax,c,0.0,0.0,0.0,false,
        meeParams,spaceCraft,thrustSunAngleConstraint,toSunVec,thrustSunAngle,panelType,eclipsing,
        false,RB,RS,desolver,reltol,abstol,maxRevs,integStep,type,returnData,writeData,
        onlyWriteDataAtSteps)
end