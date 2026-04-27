# Theory

NDBC wave records separate omnidirectional spectral density from directional
statistics. `swden` gives the scalar spectrum `S(f)`, while `swdir`, `swdir2`,
`swr1`, and `swr2` describe mean directions and Fourier spreading coefficients.
When enough directional information is available, `BuoyData.NDBC` reconstructs
unit-aware `WaveSpectra.Spectrum` objects indexed by time.

For SIRENOpt-style resource coupling, BuoyData is a data-ingestion layer. It does
not impose device assumptions. Downstream code can integrate the returned
spectrum, compute resource metrics, or feed a wave subsystem while preserving the
original station metadata and units.

Data availability and file formats are controlled by NDBC. Downloads require
network access, and requested files are cached under `~/.cache/JuliaOceanWaves` to
avoid repeated external requests during iterative simulation or optimization runs.
