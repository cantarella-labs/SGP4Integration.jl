using SGP4Integration
using SatelliteToolboxTle
using Dates


fetcher = create_tle_fetcher(CelestrakTleFetcher)
tle = fetch_tles(fetcher, satellite_name = "SATELIOT_1")[1]
to_radians(α) = α*π/180
# Kassel data:
lat = to_radians(51.299176)
lon = to_radians(9.481459)
alt = 197.777618
kassel = GeodeticPosition(lat, lon, alt)

min_elev = to_radians(15)
t0 = now(UTC)
t1 = t0 + Day(2)
passes = predict_passes(tle, kassel, t0, t1; min_elev = min_elev)
