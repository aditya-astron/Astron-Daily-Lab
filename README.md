# Astron Daily Lab

One engineering project a day across space systems, electronics and
computation: orbital mechanics, propulsion, ECLSS, comms, power, control,
signal processing and more. Each project lives in `projects/NNN-name/` with
working source, a test or verification against reference values, and a short
README.

~Astron

## Project log

| # | Date | Domain | Project | Language |
|---|------|--------|---------|----------|
| 001 | 2026-10-30 | ECLSS | [LiOH CO2 scrubber sizing for a 4-crew 30-day mission](projects/001-lioh-co2-scrubber-sizing/) | Python |
| 002 | 2026-10-31 | Telemetry framing | [CCSDS TM link: packet mux, CADU frame sync and CRC-checked demux](projects/002-ccsds-tm-frame-link/) | C |
| 003 | 2026-11-01 | Orbital mechanics | [Lambert targeting for a LEO rendezvous transfer](projects/003-lambert-rendezvous-targeting/) | Python |
| 004 | 2026-11-01 | Spacecraft power systems | [Solar array MPPT through a LEO eclipse cycle](projects/004-solar-array-mppt/) | Python |
| 005 | 2026-11-01 | Attitude determination | [Wahba's problem: TRIAD, Davenport q-method and QUEST with covariance checks](projects/005-wahba-attitude-determination/) | Rust |
| 006 | 2026-11-02 | Rocket propulsion and staging | [Optimal mass staging of a three-stage launcher (Lagrange multiplier method)](projects/006-optimal-launcher-staging/) | C++ |
| 007 | 2026-11-02 | Thermal analysis | [Transient orbital thermal model of a nadir-pointing 3U CubeSat](projects/007-cubesat-orbital-thermal/) | Python |
| 008 | 2026-11-02 | Error correction | [CCSDS K=7 R=1/2 convolutional code: Viterbi decoding with hard, 3-bit and soft decisions](projects/008-ccsds-viterbi-decoder/) | C |
| 009 | 2026-11-03 | Kalman filtering and sensor fusion | [Multiplicative EKF fusing a drifting gyro with a star tracker, checked against the Riccati steady state](projects/009-mekf-gyro-star-tracker/) | Python |
| 010 | 2026-11-03 | Laser/optical communication and link budgets | [Photon-counting PPM deep-space optical downlink: link budget, exact Poisson SER and capacity-optimal PPM order](projects/010-deep-space-optical-ppm-link/) | Python |
| 011 | 2026-11-03 | Embedded firmware | [Fixed-point Q16.16 PI speed loop for a reaction wheel with anti-windup, checked against analytic closed-loop responses](projects/011-fixed-point-wheel-speed-pi/) | C |
| 012 | 2026-11-04 | Radiation shielding | [Solar proton event dose behind aluminium vs polyethylene: Bethe stopping, CSDA transport and isotropic slab dose](projects/012-spe-proton-shielding/) | Python |
| 013 | 2026-11-04 | Navigation / Signal processing | [GPS L1 C/A Gold codes and FFT parallel code-phase acquisition with CFAR detection, checked against chi-square Pd/Pfa](projects/013-gps-ca-acquisition/) | Python |
| 014 | 2026-11-04 | ECLSS | [Oxygen loop closure: PEM electrolysis and Sabatier CO2 reduction at chemical equilibrium, with the 50% hydrogen-limited ceiling](projects/014-oxygen-loop-closure/) | Python |
| 015 | 2026-11-05 | Electronics and circuit analysis | [Damped LC input filter vs. constant-power-load instability: nodal Zout, optimal Rd-Cd damping, Middlebrook and Routh checks](projects/015-input-filter-cpl-stability/) | C++ |
| 016 | 2026-11-05 | Orbital mechanics | [Sun-synchronous repeat ground track design: J2 secular solve, sub-cycle analysis and full-cycle numerical propagation check](projects/016-sso-repeat-ground-track/) | Rust |
| 017 | 2026-11-06 | Attitude control | [Rate- and torque-limited eigenaxis slew with reaction wheels: gyroscopic cancellation, braking-curve rate loop and bang-coast-bang bound](projects/017-eigenaxis-slew-control/) | C |
| 018 | 2026-11-06 | Numerical methods | [Symplectic vs. RK4 propagation of a Molniya orbit: leapfrog and Yoshida-4 against exact Kepler f and g, with order, energy, angular momentum and reversibility checks](projects/018-symplectic-orbit-integrators/) | C++ |
| 019 | 2026-11-07 | RF link budgets | [Mars X-band downlink over a synodic cycle: DSN 34 m link budget, BI-AWGN coding threshold, VDR pass volume and conjunction outage](projects/019-mars-xband-link-volume/) | Python |
| 020 | 2026-11-07 | Spacecraft power systems | [LEO Li-ion battery sizing with a 2RC equivalent-circuit model: CC-CV charging to periodic steady state, DoD-driven string count, energy-balance and array-threshold checks](projects/020-leo-battery-sizing/) | Rust |
| 021 | 2026-11-07 | Algorithmic problem solving | [Ground-station antenna scheduling: min-cost flow on the pass time-line vs. greedy, checked against exhaustive search and the weighted-interval DP](projects/021-ground-station-pass-scheduling/) | C++ |
| 022 | 2026-11-08 | Rocket propulsion | [Minimum-length planar nozzle by the method of characteristics: wave-cancelling wall contour, mass and momentum convergence checks, Cf vs. a 15 deg wedge](projects/022-moc-minimum-length-nozzle/) | Go |
| 023 | 2026-11-08 | Error correction | [CCSDS Reed-Solomon (255,223) codec: dual-basis symbols, errors-and-erasures Berlekamp-Massey decoding, interleaving, checked against galois and the bounded-distance word error rate](projects/023-ccsds-reed-solomon-codec/) | C |
| 024 | 2026-11-09 | PID and control systems | [Coolant loop outlet-temperature PI with transport delay: SIMC vs. Ziegler-Nichols on exact-delay margins, ultimate-gain and delay-margin simulation checks, back-calculation anti-windup](projects/024-coolant-loop-pi-tuning/) | Octave |
