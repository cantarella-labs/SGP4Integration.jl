
"""
    GeodeticPosition{T}

Geodetic position with respect to the WGS84 datum.

# Fields

- `lat::T`: Latitude [rad].
- `lon::T`: Longitude [rad].
- `alt::T`: Altitude above the ellipsoid [m].
"""
struct GeodeticPosition{T}
    lat::T
    lon::T
    alt::T
end


"""
    SatellitePass{T}

A single satellite pass over a fixed ground location.

# Fields

- `aos::DateTime`: Acquisition of signal — time the satellite rises above the
  minimum elevation (UTC).
- `los::DateTime`: Loss of signal — time the satellite drops below the
  minimum elevation (UTC).
- `tca::DateTime`: Time of closest approach — time of maximum elevation (UTC).
- `max_elev::T`: Maximum elevation reached during the pass [rad].
"""
struct SatellitePass{T}
    aos::DateTime
    los::DateTime
    tca::DateTime
    max_elev::T
end
