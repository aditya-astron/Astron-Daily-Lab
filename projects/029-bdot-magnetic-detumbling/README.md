B-dot detumble, 3U CubeSat in a 500 km SSO
  inertia [0.0405, 0.0405, 0.0067] kg m^2, magnetorquers 0.2 A m^2 per axis
  orbit period 94.5 min, n = 1.1085e-3 rad/s (0.0635 deg/s)
  IGRF-2020 dipole B0 = 29806 nT, tilt 9.41 deg
  law m = -(k/|B|^2) Bdot, 1 Hz, torquers off 0.1 s around each read, 50 nT noise
  Avanzini gain k* = 2n(1+sin xi)Jmin = 2.958e-5 N m s (xi = 82.6 deg)

Reference run, w0 = (10, -8, 6) deg/s, |w0| = 14.14 deg/s
   orbits   wx deg/s   wy deg/s   wz deg/s  |w| deg/s      E (uJ)
     0.00     10.000     -8.000      6.000     14.142    1048.371
     0.25      3.992      7.031      0.201      8.087     403.253
     0.50      4.415      1.490      0.091      4.661     133.945
     0.75      2.936      0.899      0.051      3.071      58.179
     1.00      1.721      0.528      0.111      1.804      20.004
     1.50      0.710      0.125      0.136      0.733       3.221
     2.00      0.269     -0.177      0.055      0.326       0.641
     3.00     -0.075      0.127     -0.025      0.149       0.134
     4.00     -0.010      0.135     -0.043      0.143       0.116
     5.00      0.036      0.131     -0.028      0.139       0.115
     6.00      0.048      0.126     -0.011      0.135       0.112
  samples at saturation: 0.6 %
  |w| < 1.0 deg/s for good after 1.33 orbits (125 min)
  |w| < 0.5 deg/s for good after 1.73 orbits (164 min)
  residual |w| over orbit 6: mean 0.128 deg/s (2n = 0.127 deg/s)

Gain sweep (same initial state, time to |w| < 0.5 deg/s for good)
      k/k*  k/k*_max    k (N m s)     orbits
     0.250     0.041     7.396e-6       6.80
     0.500     0.083     1.479e-5       3.36
     1.000     0.165     2.958e-5       1.73
     2.000     0.331     5.917e-5       0.89
     4.000     0.662     1.183e-4       0.61
     8.000     1.323     2.367e-4       0.56
    16.000     2.647     4.733e-4       0.57
    32.000     5.294     9.467e-4       0.69
    64.000    10.588     1.893e-3       1.01

Same sweep with the ideal law (no saturation, exact rate):
      k/k*     orbits peak m (A m^2)
     0.250       5.99           0.07
     0.500       3.01           0.15
     1.000       1.53           0.28
     2.000       0.82           0.56
     4.000       0.44           1.09
     8.000       0.28           2.06
    16.000       0.68           3.73
    32.000       1.45           6.20
    64.000       2.60           9.94

Averaged theory, aligned dipole, isotropic J: slowest rate (k/J)(1 - lambda_max)
    xi deg  lambda_max      rate/(k/J)       rate at k* / n
         0      1.0000          0.0000               0.0000
        15      0.8637          0.1363               0.3431
        30      0.5917          0.4083               1.2248
        45      0.3846          0.6154               2.1012
        60      0.4327          0.5673               2.1173
        75      0.4829          0.5171               2.0332
        90      0.5000          0.5000               2.0000
```

## Reading the results

The rate does not go to zero. B-dot only drives the body to follow the
field, and in a near-polar orbit the field direction turns round about twice
per orbit, so the residual is 0.128 deg/s against 2n = 0.127 deg/s, mostly
about one axis. That is the expected hand-off state for a sun-pointing or
gravity-gradient mode, not a failure.

The ideal-law sweep shows the Avanzini effect: more gain helps only until the
decay is as fast as the field sweeps round. Past roughly k*_max (the formula
with J_max, about 8 k* for this slender body) the spin is pinned along the
field line and dragged with it, and detumbling slows again (0.28 orbits at
8 k*, 2.6 at 64 k*). With 0.2 A m^2 torquers the flight law saturates long
before that, so the optimum is flatter: anything from 4 k* to 16 k* gives
about 0.6 orbits. k* itself is a conservative choice sized on J_min that
detumbles in 1.7 orbits while saturating only 0.6 % of the time.

## Limits

Dipole field only (the full IGRF adds a few percent in direction), no
gravity-gradient, drag or residual-dipole torques, no coil inductance or
current-loop dynamics, and the magnetometer has white noise but no bias.
A calibrated bias does not matter to B-dot because it cancels in the
difference; an uncalibrated soft-iron effect from the torquers does, which is
what the quiet windows are for.
