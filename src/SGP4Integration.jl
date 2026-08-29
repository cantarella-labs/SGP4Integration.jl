module SGP4Integration
using Dates: DateTime, Period, Second
using LinearAlgebra: norm
using SatelliteToolboxPropagators: Propagators
using SatelliteToolboxTle: TLE, tle_epoch
using SatelliteToolboxTransformations:
    OrbitStateVector, PEF, TEME, geodetic_to_ecef, jd_to_date, r_eci_to_ecef

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
function predict_passes(
    tle::TLE,
    ue::GeodeticPosition,      # lat [rad], lon [rad], alt [m], WGS84
    t0::DateTime,
    t1::DateTime;  # UTC, horizon = t1 - t0
    min_elev::Float64,      # [rad]
    coarse_step::Period = Second(15),
)::Vector{SatellitePass}

    # user equipment position in ECEF frame (km)
    ue_ecef = geodetic_to_ecef(ue.lat, ue.lon, ue.alt) ./ 1e3

    #initialize the orbit propagation with SGP4 alg.
    orb = Propagators.init(Val(:SGP4), tle)

    # normal vector pointing upward on the ellipsoid.
    u = [
        cos(ue.lat)*cos(ue.lon),   # "up" at your location
        cos(ue.lat) * sin(ue.lon),
        sin(ue.lat),
    ]

    # time_offset
    t0_offset = (t0 - tle_epoch(DateTime, tle))

    # accumulation time.
    pass_duration = DateTime[]

    # elev acum and time of max elev
    in_pass = false
    aos = DateTime(0)          # only meaningful while in_pass
    max_elev = -Inf
    time_max = DateTime(0)

    # result:
    passes = SatellitePass[]
    horizon = t1 - t0
    elapsed = Second(0)

    while (elapsed < horizon)
        #propagate satellite location to next date.
        neworb = Propagators.propagate!(orb, t0_offset + elapsed, OrbitStateVector)
        new_date = jd_to_date(DateTime, neworb.t)
        # convertion matrix -> position from eci (stars reference frame) to ecef (earth center of frame)
        r_neworb_ecef = r_eci_to_ecef(TEME(), PEF(), neworb.t)

        # convert; matrix vector product (from meters to km)
        neworb_ecef = r_neworb_ecef * neworb.r ./ 1e3

        # difference between orbit vector and user equipment location
        ρ = neworb_ecef - ue_ecef
        # elevation in radians
        elev_rad = asin((ρ' * u) / norm(ρ))

        # calculate elevation
        if elev_rad > min_elev
            if !in_pass # no pass registered
                in_pass = true
                aos = new_date
                max_elev = -Inf
            end
            if elev_rad > max_elev
                max_elev = elev_rad
                time_max = new_date
            end
        elseif in_pass # end of pass
            obj = SatellitePass(aos, new_date, time_max, max_elev)
            push!(passes, obj)
            in_pass = false
        end
        elapsed += coarse_step
    end

    return passes

end
export GeodeticPosition, SatellitePass, predict_passes

end
