# 017 - Rate- and torque-limited eigenaxis slew with reaction wheels

An agile imaging satellite (J about 1200/1500/900 kg m^2 with small products of
inertia) has to retarget 90 deg using three body-axis reaction wheels rated
0.6 N m and 16 N m s, with 2.6 N m s of bias momentum already in the wheels.
The shortest rotation is about a single fixed axis (Euler's eigenaxis), but a
plain quaternion PD does not follow it: products of inertia and the gyroscopic
coupling w x (Jw + h) bend the path, and on a large error it asks the wheels
for far more torque and momentum than they have.

Goal: a feedback law that slews about the eigenaxis exactly, never exceeds a
rate set by wheel momentum or an acceleration set by wheel torque, settles to
0.01 deg, and gets close to the bang-coast-bang time bound.

## Approach

**Dynamics.** Rigid body plus wheels, no external torque:

    J w' = tau - w x (J w + h),   h' = -tau,   q' = 1/2 q (x) [w, 0]

**Limits for the chosen axis u.** Wheel i carries `w |(J u)_i|` plus at most
`|h0|`, and needs `alpha |(J u)_i|` of torque to accelerate at alpha. With an
80 % margin:

    w_max = (0.8 h_max - |h0|) / max|(J u)_i|      alpha = 0.8 tau_max / max|(J u)_i|

The leftover 20 % of torque covers the gyroscopic term `w_max |h0|` (checked).

**Control law** (after Wie, Bailey and Heiberg, JGCD 2002). With e the vector part
of the error quaternion (short way round):

    tau  = w x (J w + h) + J sat_alpha( -2k sat_L(e) - c w )
    L(e) = (c / 2k) min( sqrt(4 a |e|), w_max )

Both saturations scale whole vectors, so the commanded acceleration stays
parallel to e. Cancelling the gyroscopic term makes `w' = acc` exactly, so if
the body starts at rest, w and e stay on the initial eigenaxis the entire time.
Along that axis the law collapses to a scalar rate loop
`theta'' = -c (theta' - w_cmd)` that rides `w_cmd = min(sqrt(2 a theta), w_max)`,
the braking curve for deceleration a = 0.6 alpha; below about 9 deg it becomes
the linear PD `theta'' + c theta' + k theta = 0` with wn = 0.1 rad/s, zeta = 0.9.

**Verification** (`test_slew`):

- Hand numbers for the design: (J u)_y = 2993/sqrt(6), w_max, alpha, torque budget.
- 90 deg reference slew: axis deviation < 1e-6 rad, no torque clipping, wheel
  momentum within 80 % capacity, rate never above w_max, inertial angular
  momentum conserved to 1e-10, settled at 0.01 deg.
- The 3D angle history matches an independent RK4 integration of the scalar
  eigenaxis equation (10x finer step) to 2e-8 rad.
- A 0.02 rad slew matches the closed-form underdamped second-order response.
- Settle time is never below the bang-coast-bang bound
  `theta/w_max + w_max/alpha`, and within 30 % of it for 90 deg.
- Removing the gyroscopic term pushes the rotation axis off by about 5 deg; the
  unlimited scalar PD clips wheel torque for 155 s and demands 38 N m s.
- 40 random axes, angles (2-178 deg), bias momenta and starting attitudes all
  hold the eigenaxis, respect limits and settle.

The reference slew was also cross-checked against SciPy's DOP853 at
rtol 1e-11; the angle histories agree to 1e-6 deg.

## Build and run

    make
    ./test_slew                 # verification suite
    ./slew_demo [history.csv]   # report; optional 1 s time history

C11 and libm only.

## Sample output

```
Design (margin 80%)
  w_max = 0.4771 deg/s   alpha = 3.928e-04 rad/s^2   brake a = 2.357e-04 rad/s^2
  k = 1.0000e-02  c = 1.8000e-01  (wn = 0.100 rad/s, zeta = 0.90)
  linear PD region below 8.76 deg error

controller                settle s   axis dev    tau pk    clip s    |h| pk   rate pk
                                          deg       N m               N m s     deg/s
eigenaxis + gyro comp        258.4   1.48e-06     0.480       0.0     11.90    0.4771
eigenaxis, no gyro comp      258.9   5.13e+00     0.480       0.0     11.94    0.4766
scalar PD, unlimited         244.7   9.00e+01    13.856     155.4     37.95    1.6828

Bang-coast-bang bound with the same alpha and w_max: 209.8 s

Angle sweep about the same axis (eigenaxis law)
     deg   settle s      BCB s   ratio
       5       66.3       29.8    2.22
      30      132.6       84.1    1.58
      90      258.4      209.8    1.23
     175      436.5      388.0    1.13
```

The scalar PD only looks fast because the simulation lets its wheels hold
38 N m s; real 16 N m s wheels would saturate and lose control mid-slew. The
fixed cost of the eigenaxis law is the linear terminal phase, so its ratio to
the bound is worst on small slews and approaches 1 on long ones.
