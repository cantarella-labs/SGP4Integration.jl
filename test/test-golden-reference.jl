# test/test_golden_reference.jl
#
# Golden-reference test for SGP4Integration.predict_passes
# Reference generated with: skyfield (make_reference.py), run 2026-08-23
# Inputs frozen below — TLE, position, window, threshold must stay
# byte-identical to the values used in make_reference.py. Never re-fetch.
@testsnippet SharedData begin
#     sample_input = "Hello, World!"
#     expected_output = "Hello, World!"
# end
using Test
using Dates
using SatelliteToolboxTle

# --- frozen inputs ----------------------------------------------------------
const TLE_FROZEN = tle"""
SATELIOT_1
1 60550U 24149CL  26234.21462951  .00001807  00000+0  15683-3 0  9994
2 60550  97.6698 309.0084 0006288 198.7526 161.3466 14.98274358109951
"""

const UE       = GeodeticPosition(0.8746, 0.151531, 160.0)
const T0       = DateTime(2026, 8, 23, 12, 0, 0)        # UTC
const T1       = DateTime(2026, 8, 27, 12, 0, 0)        # UTC
const MIN_ELEV = deg2rad(21.4)
const STEP     = Second(15)

# --- expected passes (skyfield output; elevations in RADIANS — see note) ----
parse_utc(s) = DateTime(s, dateformat"yyyy-mm-ddTHH:MM:SS\Z")

const EXPECTED = [
    (parse_utc("2026-08-23T21:24:08Z"), parse_utc("2026-08-23T21:26:47Z"), parse_utc("2026-08-23T21:29:27Z"), 1.423),
    (parse_utc("2026-08-24T10:34:34Z"), parse_utc("2026-08-24T10:37:13Z"), parse_utc("2026-08-24T10:39:51Z"), 1.531),
    (parse_utc("2026-08-24T21:26:40Z"), parse_utc("2026-08-24T21:29:19Z"), parse_utc("2026-08-24T21:31:58Z"), 1.343),
    (parse_utc("2026-08-25T10:37:06Z"), parse_utc("2026-08-25T10:39:44Z"), parse_utc("2026-08-25T10:42:22Z"), 1.448),
    (parse_utc("2026-08-25T21:29:12Z"), parse_utc("2026-08-25T21:31:50Z"), parse_utc("2026-08-25T21:34:29Z"), 1.264),
    (parse_utc("2026-08-26T10:39:37Z"), parse_utc("2026-08-26T10:42:15Z"), parse_utc("2026-08-26T10:44:52Z"), 1.367),
    (parse_utc("2026-08-26T21:31:45Z"), parse_utc("2026-08-26T21:34:21Z"), parse_utc("2026-08-26T21:36:59Z"), 1.189),
    (parse_utc("2026-08-27T10:42:09Z"), parse_utc("2026-08-27T10:44:46Z"), parse_utc("2026-08-27T10:47:22Z"), 1.289),
]

# tolerances: dominated by our 15 s coarse sampling, not the physics
const TOL_EDGE = Second(30)      # AOS / LOS: 2 coarse steps
const TOL_PEAK = Second(45)      # peak time: elevation is flat near culmination
const TOL_ELEV = 0.007           # rad, ≈ 0.40°

within(a::DateTime, b::DateTime, tol::Period) = abs(a - b) <= Millisecond(tol)

end

@testitem "validation test" tags=[:unit, :validation] setup=[SharedData] begin
@testset "golden reference vs skyfield" begin
    passes = predict_passes(TLE_FROZEN, UE, T0, T1;
                            min_elev = MIN_ELEV, coarse_step = STEP)

    @test length(passes) == length(EXPECTED)

    for (p, (e_start, e_peak, e_end, e_el)) in zip(passes, EXPECTED)
        @test within(p.aos,  e_start, TOL_EDGE)
        @test within(p.los,  e_end,   TOL_EDGE)
        @test within(p.tca,  e_peak,  TOL_PEAK)
        @test abs(p.max_elev - e_el) < TOL_ELEV
    end
end

@testset "structural invariants" begin
    passes = predict_passes(TLE_FROZEN, UE, T0, T1;
                            min_elev = MIN_ELEV, coarse_step = STEP)
    for p in passes
        @test p.aos < p.tca < p.los
        @test p.max_elev >= MIN_ELEV
    end
    for i in 2:length(passes)
        @test passes[i].aos > passes[i-1].los       # ordered, non-overlapping
    end
end

@testset "empty window" begin
    @test isempty(predict_passes(TLE_FROZEN, UE,
                                 DateTime(2026, 8, 23, 12, 0, 0),
                                 DateTime(2026, 8, 23, 12, 30, 0);
                                 min_elev = MIN_ELEV, coarse_step = STEP))
end

end