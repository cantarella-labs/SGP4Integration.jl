module SGP4Integration
using Dates
using SatelliteToolboxPropagators
using SatelliteToolboxTransformations
using SatelliteToolboxTle
using LinearAlgebra




include("data_structures.jl")
"""
    predict_passes(tle, ue, t0, t1; min_elev, coarse_step = Second(15)) -> Vector{SatellitePass}

Predict the passes of the satellite described by `tle` over a ground station at `ue`
within the time window `[t0, t1]`.

The satellite is propagated with SGP4 (via `SatelliteToolboxPropagators.jl`) on a fixed
time grid of `coarse_step`. At each step the elevation of the satellite as seen from `ue`
is computed; a pass starts when the elevation rises above `min_elev` and ends when it
drops back below it.

# Arguments

- `tle::TLE`: Two-line element set of the satellite.
- `ue::GeodeticPosition`: Ground station position (lat [rad], lon [rad], alt [m], WGS84).
- `t0::DateTime`, `t1::DateTime`: Start and end of the search window (UTC).

# Keywords

- `min_elev::Float64`: Minimum elevation for the satellite to be considered visible [rad].
- `coarse_step::Period`: Sampling step of the search grid (default `Second(15)`).
  Passes shorter than one step may be missed; AOS/LOS times are only accurate to one step.

# Returns

- `Vector{SatellitePass}`: One entry per detected pass, with AOS, LOS, time of maximum
  elevation, and the maximum elevation [rad].

A pass still in progress at `t1` is discarded (no LOS inside the window).
"""
function predict_passes(tle::TLE,
               ue::GeodeticPosition,      # lat [rad], lon [rad], alt [m], WGS84
               t0::DateTime, t1::DateTime;  # UTC, horizon = t1 - t0
               min_elev::Float64,      # [rad]
               coarse_step::Period = Second(15))::Vector{SatellitePass}

    # user equipment position in ECEF frame (km)
    ue_ecef = geodetic_to_ecef(ue.lat, ue.lon, ue.alt)./1e3

    #initialize the orbit propagation with SGP4 alg.
    orb = Propagators.init(Val(:SGP4), tle)

    # normal vector pointing upward on the ellipsoid.
    u = [cos(ue.lat)*cos(ue.lon),   # "up" at your location
         cos(ue.lat)*sin(ue.lon),
         sin(ue.lat)]

    # time_offset
    t0_offset = (t0 - tle_epoch(DateTime, tle))

    # accumulation time.
    pass_duration = DateTime[]

    # elev acum and time of max elev
    max_elev = 0.0
    time_max=DateTime(1900, 1, 1)
    # result:
    passes = SatellitePass[]
    horizon = t1 - t0
    elapsed = Second(0)

    while (elapsed < horizon)
        neworb = Propagators.propagate!(orb, t0_offset + elapsed, OrbitStateVector)
        new_date = jd_to_date(DateTime, neworb.t)
        r_neworb_ecef = r_eci_to_ecef(TEME(), PEF(), date_to_jd(new_date))
        neworb_ecef = r_neworb_ecef * neworb.r ./ 1e3
        ρ = neworb_ecef - ue_ecef
        elev_rad = asin((ρ' * u) / norm(ρ))
        if elev_rad > min_elev
            if length(pass_duration) < 1 # no pass registered
                push!(pass_duration, new_date)
            elseif length(pass_duration) == 1 #during pass
                if elev_rad > max_elev
                    max_elev = elev_rad
                    time_max = new_date
                end
            end
        else
            if length(pass_duration) == 1 # end of pass
                push!(pass_duration, new_date)
                obj = SatellitePass(
                    pass_duration[1],
                    pass_duration[2],
                    time_max,
                    max_elev
                )
                push!(passes, obj)
                pass_duration = DateTime[]
                max_elev = 0.0
                time_max = DateTime(1900, 1, 1)
            end
        end
        elapsed += coarse_step
    end

    return passes
    
end
export GeodeticPosition, SatellitePass, predict_passes

end