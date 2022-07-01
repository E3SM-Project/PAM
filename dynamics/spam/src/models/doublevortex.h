
struct DoubleVortex {
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 5000000._fp;
  static real constexpr Ly = 5000000._fp;
  static real constexpr coriolis = 0.00006147_fp;
  static real constexpr H0 = 750.0_fp;
  static real constexpr ox = 0.1_fp;
  static real constexpr oy = 0.1_fp;
  static real constexpr sigmax = 3._fp / 40._fp * Lx;
  static real constexpr sigmay = 3._fp / 40._fp * Ly;
  static real constexpr dh = 75.0_fp;
  static real constexpr xc1 = (0.5_fp - ox) * Lx;
  static real constexpr yc1 = (0.5_fp - oy) * Ly;
  static real constexpr xc2 = (0.5_fp + ox) * Lx;
  static real constexpr yc2 = (0.5_fp + oy) * Ly;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr yc = 0.5_fp * Ly;
  static real constexpr c = 0.05_fp;
  static real constexpr a = 1.0_fp / 3.0_fp;
  static real constexpr D = 0.5_fp * Lx;

  static real YAKL_INLINE coriolis_f(real x, real y) { return coriolis; }

  static real YAKL_INLINE h_f(real x, real y) {
    real xprime1 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc1));
    real yprime1 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc1));
    real xprime2 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc2));
    real yprime2 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc2));
    real xprimeprime1 =
        Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc1));
    real yprimeprime1 =
        Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc1));
    real xprimeprime2 =
        Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc2));
    real yprimeprime2 =
        Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc2));

    return H0 - dh * (exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1)) +
                      exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2)) -
                      4._fp * pi * sigmax * sigmay / Lx / Ly);
  }

  static vec<2> YAKL_INLINE v_f(real x, real y) {
    vec<2> vvec;

    real xprime1 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc1));
    real yprime1 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc1));
    real xprime2 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc2));
    real yprime2 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc2));
    real xprimeprime1 =
        Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc1));
    real yprimeprime1 =
        Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc1));
    real xprimeprime2 =
        Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc2));
    real yprimeprime2 =
        Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc2));

    vvec.u =
        -g * dh / coriolis / sigmay *
        (yprimeprime1 * exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1)) +
         yprimeprime2 * exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2)));
    vvec.v =
        g * dh / coriolis / sigmax *
        (xprimeprime1 * exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1)) +
         xprimeprime2 * exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2)));
    return vvec;
  }

  static real YAKL_INLINE S_f(real x, real y) {
    // real sval = g * (1. + c * sin(2. * M_PI / Lx * (x - xc)) * sin(2. * M_PI
    // / Ly * (y - yc)) * exp(-((x-xc)*(x-xc) + (y-yc)*(y-yc))/(a*a*D*D)));
    real sval =
        g * (1._fp + c * exp(-((x - xc) * (x - xc) + (y - yc) * (y - yc)) /
                             (a * a * D * D)));
    // real sval = g * (1. + c * sin(2. * M_PI / Lx * (x- xc)));
    // real sval = g;
    // real sval = g * (1. + c * ((x > 0.35 * Lx && x < 0.65 * Lx && y > 0.35 *
    // Ly
    // && y < 0.65 * Ly ) ? 1. : 0.));
    return sval * h_f(x, y);
  }
};
