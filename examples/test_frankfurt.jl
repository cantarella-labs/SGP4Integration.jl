using Dates
using SatelliteToolboxPropagators
using SatelliteToolboxTransformations
using SatelliteToolboxTle


fetcher = create_tle_fetcher(CelestrakTleFetcher)
tle = fetch_tles(fetcher, satellite_name = "SATELIOT_1")[1]
tle_epoch(DateTime, tle)
orb = Propagators.init(Val(:SGP4), tle)
neworb = Propagators.propagate!(orb, Dates.Minute(60), OrbitStateVector)
new_date = jd_to_date(DateTime, neworb.t)


"Holds the geodetic position with regards to the wgs84 datum."
struct GeodeticPosition{T}
    "latitude in radians"
    lat::T
    "longitude in radians"
    lon::T
    "altitude in m"
    alt::T
end

frank = GeodeticPosition(0.8746, 0.151531, 160.0)

frank_ecef = geodetic_to_ecef(frank.lat, frank.lon, frank.alt) ./ 1e3
r_neworb_ecef = r_eci_to_ecef(TEME(), PEF(), date_to_jd(new_date))
neworb_ecef = r_neworb_ecef * neworb.r ./ 1e3

using LinearAlgebra

u_ = frank_ecef ./ norm(frank_ecef)
ρ = neworb_ecef - frank_ecef
u = [
    cos(frank.lat)*cos(frank.lon),   # "up" at your location
    cos(frank.lat) * sin(frank.lon),
    sin(frank.lat),
]
elev_rad = asin((ρ' * u) / norm(ρ))

## Pass prediction alg. → very stupid.
# calculate the orbit for the next 4 days and check when the next time the sat passes next to you.

days = 4
days_in_secods = days*24*3600
acc_time = 15
step = 15
min_elev = (21.4*π)/180
pass_duration = DateTime[]
passes = Vector{DateTime}[]
while (acc_time < days_in_secods)
    neworb = Propagators.propagate!(orb, Dates.Second(acc_time), OrbitStateVector)
    new_date = jd_to_date(DateTime, neworb.t)
    r_neworb_ecef = r_eci_to_ecef(TEME(), PEF(), date_to_jd(new_date))
    neworb_ecef = r_neworb_ecef * neworb.r ./ 1e3
    ρ = neworb_ecef - frank_ecef
    elev_rad = asin((ρ' * u) / norm(ρ))
    if elev_rad > min_elev
        if length(pass_duration) < 1 # no pass registered
            push!(pass_duration, new_date)
        end
    else
        if length(pass_duration) == 1 # end of pass
            push!(pass_duration, new_date)
            push!(passes, copy(pass_duration))
            pass_duration = DateTime[]
        end
    end
    acc_time += step
end
