/*
 * Astron Built this AP
 * Project 039 | Rocket propulsion | BATES solid motor internal ballistics
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-15 IST
 *
 * Sizes the throat of a four-segment BATES motor for a 5 MPa peak chamber
 * pressure, simulates ignition fill, burn and tail-off, and prints the
 * performance summary. `./motor --csv` writes the pressure/thrust trace.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ballistics.h"
#include "motor_design.h"

#define TRACE_CAP 4096

int main(int argc, char **argv)
{
    int csv = argc > 1 && strcmp(argv[1], "--csv") == 0;
    Propellant p = design_propellant();
    Grain g = design_grain();

    Nozzle nz;
    nz.At  = size_throat(&p, &g, P_DESIGN_MAX);
    /* Expand to ambient at the average chamber pressure, ~0.95 of peak for
     * a near-neutral grain; this keeps the burn free of overexpansion loss
     * worth worrying about at sea level. */
    nz.eps = eps_for_optimum(0.95 * P_DESIGN_MAX, PA_SEA_LEVEL, p.gamma);

    static TraceSample tr[TRACE_CAP];
    MotorResult r = simulate_motor(&p, &g, &nz, PA_SEA_LEVEL, 1e-5,
                                   P_IGNITER, tr, TRACE_CAP, 2000);

    if (csv) {
        printf("t_s,web_mm,Pc_MPa,F_N,Ab_cm2\n");
        for (size_t i = 0; i < r.n_samples; ++i)
            printf("%.4f,%.4f,%.5f,%.2f,%.3f\n", tr[i].t, tr[i].x * 1e3,
                   tr[i].Pc / 1e6, tr[i].F, tr[i].Ab * 1e4);
        return 0;
    }

    double w = grain_web(&g);
    double A0 = grain_burn_area(&g, 0.0), Aw = grain_burn_area(&g, w * (1 - 1e-12));
    double xs = (g.L - 2.0 * g.d) / 6.0, Am = grain_burn_area(&g, xs);
    double Dt = sqrt(4.0 * nz.At / M_PI);
    double cfs = thrust_coeff(nz.eps, p.gamma, P_DESIGN_MAX, PA_SEA_LEVEL);
    double cfv = thrust_coeff(nz.eps, p.gamma, 1.0, 0.0);
    double g0 = 9.80665;

    printf("BATES motor, %d x (D %.1f mm, d %.1f mm, L %.1f mm)\n",
           g.N, g.D * 1e3, g.d * 1e3, g.L * 1e3);
    printf("  web %.2f mm, propellant %.3f kg, volumetric loading %.3f\n",
           w * 1e3, r.m_prop,
           grain_prop_volume(&g, 0) / (0.25 * M_PI * g.D * g.D * g.Lcase));
    printf("  Ab start/peak/end  %.1f / %.1f / %.1f cm^2  (peak at x = %.2f mm,"
           " +%.2f %%)\n", A0 * 1e4, Am * 1e4, Aw * 1e4, xs * 1e3,
           100.0 * (Am / A0 - 1.0));
    printf("Gas: c* %.1f m/s (ideal %.1f), r(5 MPa) %.2f mm/s\n",
           cstar(&p), cstar_ideal(&p), p.a * pow(5e6, p.n) * 1e3);
    printf("Nozzle: throat %.2f mm dia, Kn %.0f..%.0f, eps %.2f, "
           "pe/pc %.4f\n", Dt * 1e3, A0 / nz.At, Am / nz.At, nz.eps,
           pe_over_pc(nz.eps, p.gamma));
    printf("  CF sea level at MEOP %.4f, CF vacuum %.4f\n", cfs, cfv);
    printf("Simulation (dt 10 us, ignition at %.1f MPa)\n", P_IGNITER / 1e6);
    printf("  peak Pc %.4f MPa (target %.1f), peak thrust %.0f N\n",
           r.P_max / 1e6, P_DESIGN_MAX / 1e6, r.F_max);
    printf("  web burnout %.3f s at %.3f MPa, run ends %.3f s (Pc < 2 pa)\n",
           r.t_burnout, r.P_burnout / 1e6, r.t_end);
    printf("  total impulse  sea level %.0f N s, vacuum %.0f N s\n",
           r.I_total, r.I_vac);
    printf("  delivered Isp  sea level %.1f s, vacuum %.1f s\n",
           r.I_total / (r.m_prop * g0), r.I_vac / (r.m_prop * g0));
    printf("  mass balance   nozzle %.5f + chamber %.5f = %.5f kg vs %.5f kg"
           " loaded\n", r.m_nozzle, r.m_gas_end, r.m_nozzle + r.m_gas_end,
           r.m_prop);
    double I_class = r.I_total;
    const char *cls = "?";
    const char *names = "ABCDEFGHIJKLMNO";
    for (int i = 0; i < 15; ++i)
        if (I_class <= 2.5 * pow(2.0, i)) { cls = &names[i]; break; }
    printf("  impulse class  %c (avg thrust %.0f N)\n", *cls,
           r.I_total / r.t_end);
    return 0;
}
