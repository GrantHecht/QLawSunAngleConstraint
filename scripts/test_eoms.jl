using DrWatson
@quickactivate "QLawSunAngleConstraint"

# Here you may include files from the source directory
include(srcdir("include.jl"))
using DataFrames, CSV
import Interpolations

function main()
    # Define initial state and epoch
    furnsh_kernals()
    init_epoch = utc2et("2000-03-22T00:00:00")
    mee0       = SVector(11359.07, 0.7306, 0.0, 0.2539676, 0.0, 0.0)

    # Create the serialized ephemeris
    UTC0        = "2000-03-18T00:00:00"
    UTCF        = "2020-04-10T00:00:00"
    n_points    = 10000
    targ_ids    = [10, 301,]
    obs_id      = 399
    reff        = "J2000"
    ephem = create_serialized_ephemerides(
        UTC0, UTCF, n_points, targ_ids, obs_id, reff,
    )
    unload_kernals()

    # Construct ForceModel
    central_body = 399
    fm = ForceModel{Ephemerides}(
        init_epoch,
        get_gm(central_body),
        true,
        true,
        ephem,
    )

    # Propagate dynamics
    fun  = ODEFunction{false}((u,p,t) -> mee_equations_of_motion(u,p,t,fm))
    prob = ODEProblem(fun, mee0, (0.0, 10*86400.0))
    sol =  solve(prob, Tsit5())

    # Load GMAT ephemeris for same propagation
    gmat_ephem = DataFrame(
        CSV.File(
            joinpath(@__DIR__, "..", "data", "gmat_ephem_data", "report_sun_luna.txt");
            delim = ' ',
            ignorerepeated = true,
        )
    )


    # Plot results
    f = Figure(resolution = (1000, 1200))
    ax_p = Axis(f[1, 1], xlabel = "Time, sec", ylabel = "Semi-perimiter, km")
    pdiff = begin
        gmatp = Interpolations.interpolate(
            (gmat_ephem[:, 1],), 
            gmat_ephem[:, 2], 
            Interpolations.Gridded(Interpolations.Linear()),
        )
        [sol[i][1] - gmatp(sol.t[i]) for i in eachindex(sol)]
    end
    lines!(ax_p, sol.t, pdiff)
    ax_f = Axis(f[2, 1], xlabel = "Time, sec", ylabel = "f, n.d.")
    fdiff = begin
        gmatf = Interpolations.interpolate(
            (gmat_ephem[:, 1],), 
            gmat_ephem[:, 3], 
            Interpolations.Gridded(Interpolations.Linear()),
        )
        [sol[i][2] - gmatf(sol.t[i]) for i in eachindex(sol)]
    end
    lines!(ax_f, sol.t, fdiff)
    ax_g = Axis(f[3, 1], xlabel = "Time, sec", ylabel = "g, n.d.")
    gdiff = begin
        gmatg = Interpolations.interpolate(
            (gmat_ephem[:, 1],), 
            gmat_ephem[:, 4], 
            Interpolations.Gridded(Interpolations.Linear()),
        )
        [sol[i][3] - gmatg(sol.t[i]) for i in eachindex(sol)]
    end
    lines!(ax_g, sol.t, gdiff)
    ax_h = Axis(f[4, 1], xlabel = "Time, sec", ylabel = "h, n.d.")
    hdiff = begin
        gmath = Interpolations.interpolate(
            (gmat_ephem[:, 1],), 
            gmat_ephem[:, 5], 
            Interpolations.Gridded(Interpolations.Linear()),
        )
        [sol[i][4] - gmath(sol.t[i]) for i in eachindex(sol)]
    end
    lines!(ax_h, sol.t, hdiff)
    ax_k = Axis(f[5, 1], xlabel = "Time, sec", ylabel = "k, n.d.")
    kdiff = begin
        gmatk = Interpolations.interpolate(
            (gmat_ephem[:, 1],), 
            gmat_ephem[:, 6], 
            Interpolations.Gridded(Interpolations.Linear()),
        )
        [sol[i][5] - gmatk(sol.t[i]) for i in eachindex(sol)]
    end
    lines!(ax_k, sol.t, kdiff)

    return f
end

main()