/*
 * Astron Built this AP
 * Project 011 | Embedded firmware | Fixed-point PI speed loop for a reaction wheel with anti-windup
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-03 IST
 */
#include "wheel.h"

#include <math.h>

/* w(t) = w_inf + (w0 - w_inf) exp(-b t / J),  w_inf = (Kt i - tau) / b */
void wheel_step(wheel_t *wh, double i, double tau_load, double dt)
{
    double w_inf = (wh->kt * i - tau_load) / wh->b;
    wh->w = w_inf + (wh->w - w_inf) * exp(-wh->b * dt / wh->J);
}
