# ==== This script can be used to generate serialized ephemerides for the EOMs

using DrWatson
@quickactivate "QLawSunAngleConstraint"

# Here you may include files from the source directory
include(srcdir("include.jl"))

# Furnsh kernels
furnsh_kernals()

# Create the serialized ephemeris
UTC0        = "2000-03-22T00:00:00"
UTCF        = "2020-04-10T00:00:00"
n_points    = 10000
targ_ids    = [10, 301,]
obs_id      = 399
reff        = "J2000"
create_serialized_ephemerides(
    UTC0, UTCF, n_points, targ_ids, obs_id, reff,
)
