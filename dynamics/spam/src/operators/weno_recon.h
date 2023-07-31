#pragma once

#include "common.h"

namespace pamc {
real YAKL_INLINE interp_weno3(real phim1, real phi, real phip1) {
  const real p0 = (-1.0_fp / 2.0_fp) * phim1 + (3.0_fp / 2.0_fp) * phi;
  const real p1 = (1.0_fp / 2.0_fp) * phi + (1.0_fp / 2.0_fp) * phip1;

  const real beta1 = (phip1 - phi) * (phip1 - phi);
  const real beta0 = (phi - phim1) * (phi - phim1);

  const real alpha0 =
      (1.0_fp / 3.0_fp) / ((beta0 + 1e-10_fp) * (beta0 + 1.0e-10_fp));
  const real alpha1 =
      (2.0_fp / 3.0_fp) / ((beta1 + 1e-10_fp) * (beta1 + 1.0e-10_fp));

  const real alpha_sum_inv = 1.0_fp / (alpha0 + alpha1);

  const real w0 = alpha0 * alpha_sum_inv;
  const real w1 = alpha1 * alpha_sum_inv;

  return w0 * p0 + w1 * p1;
};

real YAKL_INLINE interp_weno5(real phim2, real phim1, real phi, real phip1,
                              real phip2) {

  const real p0 = (1.0_fp / 3.0_fp) * phim2 - (7.0_fp / 6.0_fp) * phim1 +
                  (11.0_fp / 6.0_fp) * phi;
  const real p1 = (-1.0_fp / 6.0_fp) * phim1 + (5.0_fp / 6.0_fp) * phi +
                  (1.0_fp / 3.0_fp) * phip1;
  const real p2 = (1.0_fp / 3.0_fp) * phi + (5.0_fp / 6.0_fp) * phip1 -
                  (1.0_fp / 6.0_fp) * phip2;

  const real beta2 = (13.0_fp / 12.0_fp * (phi - 2.0_fp * phip1 + phip2) *
                          (phi - 2.0_fp * phip1 + phip2) +
                      0.25_fp * (3.0_fp * phi - 4.0_fp * phip1 + phip2) *
                          (3.0_fp * phi - 4.0_fp * phip1 + phip2));
  const real beta1 = (13.0_fp / 12.0_fp * (phim1 - 2.0_fp * phi + phip1) *
                          (phim1 - 2.0_fp * phi + phip1) +
                      0.25_fp * (phim1 - phip1) * (phim1 - phip1));
  const real beta0 = (13.0_fp / 12.0_fp * (phim2 - 2.0_fp * phim1 + phi) *
                          (phim2 - 2.0_fp * phim1 + phi) +
                      0.25_fp * (phim2 - 4.0_fp * phim1 + 3.0_fp * phi) *
                          (phim2 - 4.0_fp * phim1 + 3.0_fp * phi));

  const real alpha0 = 0.1_fp / ((beta0 + 1e-10_fp) * (beta0 + 1e-10_fp));
  const real alpha1 = 0.6_fp / ((beta1 + 1e-10_fp) * (beta1 + 1e-10_fp));
  const real alpha2 = 0.3_fp / ((beta2 + 1e-10_fp) * (beta2 + 1e-10_fp));

  const real alpha_sum_inv = 1.0_fp / (alpha0 + alpha1 + alpha2);

  const real w0 = alpha0 * alpha_sum_inv;
  const real w1 = alpha1 * alpha_sum_inv;
  const real w2 = alpha2 * alpha_sum_inv;

  return w0 * p0 + w1 * p1 + w2 * p2;
};

real YAKL_INLINE interp_weno7(real phim3, real phim2, real phim1, real phi,
                              real phip1, real phip2, real phip3) {

  const real p0 = (-1.0_fp / 4.0_fp) * phim3 + (13.0_fp / 12.0_fp) * phim2 +
                  (-23.0_fp / 12.0_fp) * phim1 + (25.0_fp / 12.0_fp) * phi;
  const real p1 = (1.0_fp / 12.0_fp) * phim2 + (-5.0_fp / 12.0_fp) * phim1 +
                  (13.0_fp / 12.0_fp) * phi + (1.0_fp / 4.0_fp) * phip1;
  const real p2 = (-1.0_fp / 12.0_fp) * phim1 + (7.0_fp / 12.0_fp) * phi +
                  (7.0_fp / 12.0_fp) * phip1 + (-1.0_fp / 12.0_fp) * phip2;
  const real p3 = (1.0_fp / 4.0_fp) * phi + (13.0_fp / 12.0_fp) * phip1 +
                  (-5.0_fp / 12.0_fp) * phip2 + (1.0_fp / 12.0_fp) * phip3;

  const real beta0 =
      (phim3 * (547.0_fp * phim3 - 3882.0_fp * phim2 + 4642.0_fp * phim1 -
                1854.0_fp * phi) +
       phim2 * (7043.0_fp * phim2 - 17246.0_fp * phim1 + 7042.0_fp * phi) +
       phim1 * (11003.0_fp * phim1 - 9402.0_fp * phi) + 2107.0_fp * phi * phi);
  const real beta1 =
      (phim2 * (267.0_fp * phim2 - 1642.0_fp * phim1 + 1602.0_fp * phi -
                494.0_fp * phip1) +
       phim1 * (2843.0_fp * phim1 - 5966.0_fp * phi + 1922.0_fp * phip1) +
       phi * (3443.0_fp * phi - 2522.0_fp * phip1) + 547.0_fp * phip1 * phip1);
  const real beta2 =
      (phim1 * (547.0_fp * phim1 - 2522.0_fp * phi + 1922.0_fp * phip1 -
                494.0_fp * phip2) +
       phi * (3443.0_fp * phi - 5966.0_fp * phip1 + 1602.0_fp * phip2) +
       phip1 * (2843.0_fp * phip1 - 1642.0_fp * phip2) +
       267.0_fp * phip2 * phip2);
  const real beta3 =
      (phi * (2107.0_fp * phi - 9402.0_fp * phip1 + 7042.0_fp * phip2 -
              1854.0_fp * phip3) +
       phip1 * (11003.0_fp * phip1 - 17246.0_fp * phip2 + 4642.0_fp * phip3) +
       phip2 * (7043.0_fp * phip2 - 3882.0_fp * phip3) +
       547.0_fp * phip3 * phip3);

  const real alpha0 =
      (1.0_fp / 35.0_fp) / ((beta0 + 1e-10_fp) * (beta0 + 1e-10_fp));
  const real alpha1 =
      (12.0_fp / 35.0_fp) / ((beta1 + 1e-10_fp) * (beta1 + 1e-10_fp));
  const real alpha2 =
      (18.0_fp / 35.0_fp) / ((beta2 + 1e-10_fp) * (beta2 + 1e-10_fp));
  const real alpha3 =
      (4.0_fp / 35.0_fp) / ((beta3 + 1e-10_fp) * (beta3 + 1e-10_fp));

  const real alpha_sum_inv = 1.0_fp / (alpha0 + alpha1 + alpha2 + alpha3);

  const real w0 = alpha0 * alpha_sum_inv;
  const real w1 = alpha1 * alpha_sum_inv;
  const real w2 = alpha2 * alpha_sum_inv;
  const real w3 = alpha3 * alpha_sum_inv;

  return w0 * p0 + w1 * p1 + w2 * p2 + w3 * p3;
};

real YAKL_INLINE interp_weno9(real phim4, real phim3, real phim2, real phim1,
                              real phi, real phip1, real phip2, real phip3,
                              real phip4) {
  const real p0 = (1.0_fp / 5.0_fp) * phim4 + (-21.0_fp / 20.0_fp) * phim3 +
                  (137.0_fp / 60.0_fp) * phim2 + (-163.0_fp / 60.0_fp) * phim1 +
                  (137.0_fp / 60.0_fp) * phi;
  const real p1 = (-1.0_fp / 20.0_fp) * phim3 + (17.0_fp / 60.0_fp) * phim2 +
                  (-43.0_fp / 60.0_fp) * phim1 + (77.0_fp / 60.0_fp) * phi +
                  (1.0_fp / 5.0_fp) * phip1;
  const real p2 = (1.0_fp / 30.0_fp) * phim2 + (-13.0_fp / 60.0_fp) * phim1 +
                  (47.0_fp / 60.0_fp) * phi + (9.0_fp / 20.0_fp) * phip1 +
                  (-1.0_fp / 20.0_fp) * phip2;
  const real p3 = (-1.0_fp / 20.0_fp) * phim1 + (9.0_fp / 20.0_fp) * phi +
                  (47.0_fp / 60.0_fp) * phip1 + (-13.0_fp / 60.0_fp) * phip2 +
                  (1.0_fp / 30.0_fp) * phip3;
  const real p4 = (1.0_fp / 5.0_fp) * phi + (77.0_fp / 60.0_fp) * phip1 +
                  (-43.0_fp / 60.0_fp) * phip2 + (17.0_fp / 60.0_fp) * phip3 +
                  (-1.0_fp / 20.0_fp) * phip4;

  const real beta0 =
      (phim4 * (22658.0_fp * phim4 - 208501.0_fp * phim3 + 364863.0_fp * phim2 -
                288007.0_fp * phim1 + 86329.0_fp * phi) +
       phim3 * (482963.0_fp * phim3 - 1704396.0_fp * phim2 +
                1358458.0_fp * phim1 - 411487.0_fp * phi) +
       phim2 *
           (1521393.0_fp * phim2 - 2462076.0_fp * phim1 + 758823.0_fp * phi) +
       phim1 * (1020563.0_fp * phim1 - 649501.0_fp * phi) +
       107918.0_fp * phi * phi);
  const real beta1 =
      (phim3 * (6908.0_fp * phim3 - 60871.0_fp * phim2 + 99213.0_fp * phim1 -
                70237.0_fp * phi + 18079.0_fp * phip1) +
       phim2 * (138563.0_fp * phim2 - 464976.0_fp * phim1 + 337018.0_fp * phi -
                88297.0_fp * phip1) +
       phim1 * (406293.0_fp * phim1 - 611976.0_fp * phi + 165153.0_fp * phip1) +
       phi * (242723.0_fp * phi - 140251.0_fp * phip1) +
       22658.0_fp * phip1 * phip1);
  const real beta2 =
      (phim2 * (6908.0_fp * phim2 - 51001.0_fp * phim1 + 67923.0_fp * phi -
                38947.0_fp * phip1 + 8209.0_fp * phip2) +
       phim1 * (104963.0_fp * phim1 - 299076.0_fp * phi + 179098.0_fp * phip1 -
                38947.0_fp * phip2) +
       phi * (231153.0_fp * phi - 299076.0_fp * phip1 + 67923.0_fp * phip2) +
       phip1 * (104963.0_fp * phip1 - 51001.0_fp * phip2) +
       6908.0_fp * phip2 * phip2);
  const real beta3 =
      (phim1 * (22658.0_fp * phim1 - 140251.0_fp * phi + 165153.0_fp * phip1 -
                88297.0_fp * phip2 + 18079.0_fp * phip3) +
       phi * (242723.0_fp * phi - 611976.0_fp * phip1 + 337018.0_fp * phip2 -
              70237.0_fp * phip3) +
       phip1 *
           (406293.0_fp * phip1 - 464976.0_fp * phip2 + 99213.0_fp * phip3) +
       phip2 * (138563.0_fp * phip2 - 60871.0_fp * phip3) +
       6908.0_fp * phip3 * phip3);
  const real beta4 =
      (phi * (107918.0_fp * phi - 649501.0_fp * phip1 + 758823.0_fp * phip2 -
              411487.0_fp * phip3 + 86329.0_fp * phip4) +
       phip1 * (1020563.0_fp * phip1 - 2462076.0_fp * phip2 +
                1358458.0_fp * phip3 - 288007.0_fp * phip4) +
       phip2 *
           (1521393.0_fp * phip2 - 1704396.0_fp * phip3 + 364863.0_fp * phip4) +
       phip3 * (482963.0_fp * phip3 - 208501.0_fp * phip4) +
       22658.0_fp * phip4 * phip4);

  const real alpha0 =
      (1.0_fp / 126.0_fp) / ((beta0 + 1e-10_fp) * (beta0 + 1e-10_fp));
  const real alpha1 =
      (10.0_fp / 63.0_fp) / ((beta1 + 1e-10_fp) * (beta1 + 1e-10_fp));
  const real alpha2 =
      (10.0_fp / 21.0_fp) / ((beta2 + 1e-10_fp) * (beta2 + 1e-10_fp));
  const real alpha3 =
      (20.0_fp / 63.0_fp) / ((beta3 + 1e-10_fp) * (beta3 + 1e-10_fp));
  const real alpha4 =
      (5.0_fp / 126.0_fp) / ((beta4 + 1e-10_fp) * (beta4 + 1e-10_fp));

  const real alpha_sum_inv =
      1.0_fp / (alpha0 + alpha1 + alpha2 + alpha3 + alpha4);

  const real w0 = alpha0 * alpha_sum_inv;
  const real w1 = alpha1 * alpha_sum_inv;
  const real w2 = alpha2 * alpha_sum_inv;
  const real w3 = alpha3 * alpha_sum_inv;
  const real w4 = alpha4 * alpha_sum_inv;

  return w0 * p0 + w1 * p1 + w2 * p2 + w3 * p3 + w4 * p4;
};

real YAKL_INLINE interp_weno11(real phim5, real phim4, real phim3, real phim2,
                               real phim1, real phi, real phip1, real phip2,
                               real phip3, real phip4, real phip5) {

  const real p0 = ((-1.0_fp / 6.0_fp) * phim5 + (31.0_fp / 30.0_fp) * phim4 +
                   (-163.0_fp / 60.0_fp) * phim3 + (79.0_fp / 20.0_fp) * phim2 +
                   (-71.0_fp / 20.0_fp) * phim1 + (49.0_fp / 20.0_fp) * phi);
  const real p1 = ((1.0_fp / 30.0_fp) * phim4 + (-13.0_fp / 60.0_fp) * phim3 +
                   (37.0_fp / 60.0_fp) * phim2 + (-21.0_fp / 20.0_fp) * phim1 +
                   (29.0_fp / 20.0_fp) * phi + (1.0_fp / 6.0_fp) * phip1);
  const real p2 = ((-1.0_fp / 60.0_fp) * phim3 + (7.0_fp / 60.0_fp) * phim2 +
                   (-23.0_fp / 60.0_fp) * phim1 + (19.0_fp / 20.0_fp) * phi +
                   (11.0_fp / 30.0_fp) * phip1 + (-1.0_fp / 30.0_fp) * phip2);
  const real p3 = ((1.0_fp / 60.0_fp) * phim2 + (-2.0_fp / 15.0_fp) * phim1 +
                   (37.0_fp / 60.0_fp) * phi + (37.0_fp / 60.0_fp) * phip1 +
                   (-2.0_fp / 15.0_fp) * phip2 + (1.0_fp / 60.0_fp) * phip3);
  const real p4 = ((-1.0_fp / 30.0_fp) * phim1 + (11.0_fp / 30.0_fp) * phi +
                   (19.0_fp / 20.0_fp) * phip1 + (-23.0_fp / 60.0_fp) * phip2 +
                   (7.0_fp / 60.0_fp) * phip3 + (-1.0_fp / 60.0_fp) * phip4);
  const real p5 = ((1.0_fp / 6.0_fp) * phi + (29.0_fp / 20.0_fp) * phip1 +
                   (-21.0_fp / 20.0_fp) * phip2 + (37.0_fp / 60.0_fp) * phip3 +
                   (-13.0_fp / 60.0_fp) * phip4 + (1.0_fp / 30.0_fp) * phip5);

  const real beta0 = (phim5 * (1152561.0_fp * phim5 - 12950184.0_fp * phim4 +
                               29442256.0_fp * phim3 - 33918804.0_fp * phim2 +
                               19834350.0_fp * phim1 - 4712740.0_fp * phi) +
                      phim4 * (36480687.0_fp * phim4 - 166461044.0_fp * phim3 +
                               192596472.0_fp * phim2 - 113206788.0_fp * phim1 +
                               27060170.0_fp * phi) +
                      phim3 * (190757572.0_fp * phim3 - 444003904.0_fp * phim2 +
                               262901672.0_fp * phim1 - 63394124.0_fp * phi) +
                      phim2 * (260445372.0_fp * phim2 - 311771244.0_fp * phim1 +
                               76206736.0_fp * phi) +
                      phim1 * (94851237.0_fp * phim1 - 47460464.0_fp * phi) +
                      6150211.0_fp * phi * phi);

  const real beta1 = (phim4 * (271779.0_fp * phim4 - 3015728.0_fp * phim3 +
                               6694608.0_fp * phim2 - 7408908.0_fp * phim1 +
                               4067018.0_fp * phi - 880548.0_fp * phip1) +
                      phim3 * (8449957.0_fp * phim3 - 37913324.0_fp * phim2 +
                               42405032.0_fp * phim1 - 23510468.0_fp * phi +
                               5134574.0_fp * phip1) +
                      phim2 * (43093692.0_fp * phim2 - 97838784.0_fp * phim1 +
                               55053752.0_fp * phi - 12183636.0_fp * phip1) +
                      phim1 * (56662212.0_fp * phim1 - 65224244.0_fp * phi +
                               14742480.0_fp * phip1) +
                      phi * (19365967.0_fp * phi - 9117992.0_fp * phip1) +
                      1152561.0_fp * phip1 * phip1);

  const real beta2 = (phim3 * (139633.0_fp * phim3 - 1429976.0_fp * phim2 +
                               2863984.0_fp * phim1 - 2792660.0_fp * phi +
                               1325006.0_fp * phip1 - 245620.0_fp * phip2) +
                      phim2 * (3824847.0_fp * phim2 - 15880404.0_fp * phim1 +
                               15929912.0_fp * phi - 7727988.0_fp * phip1 +
                               1458762.0_fp * phip2) +
                      phim1 * (17195652.0_fp * phim1 - 35817664.0_fp * phi +
                               17905032.0_fp * phip1 - 3462252.0_fp * phip2) +
                      phi * (19510972.0_fp * phi - 20427884.0_fp * phip1 +
                             4086352.0_fp * phip2) +
                      phip1 * (5653317.0_fp * phip1 - 2380800.0_fp * phip2) +
                      271779.0_fp * phip2 * phip2);

  const real beta3 = (phim2 * (271779.0_fp * phim2 - 2380800.0_fp * phim1 +
                               4086352.0_fp * phi - 3462252.0_fp * phip1 +
                               1458762.0_fp * phip2 - 245620.0_fp * phip3) +
                      phim1 * (5653317.0_fp * phim1 - 20427884.0_fp * phi +
                               17905032.0_fp * phip1 - 7727988.0_fp * phip2 +
                               1325006.0_fp * phip3) +
                      phi * (19510972.0_fp * phi - 35817664.0_fp * phip1 +
                             15929912.0_fp * phip2 - 2792660.0_fp * phip3) +
                      phip1 * (17195652.0_fp * phip1 - 15880404.0_fp * phip2 +
                               2863984.0_fp * phip3) +
                      phip2 * (3824847.0_fp * phip2 - 1429976.0_fp * phip3) +
                      139633.0_fp * phip3 * phip3);

  const real beta4 = (phim1 * (1152561.0_fp * phim1 - 9117992.0_fp * phi +
                               14742480.0_fp * phip1 - 12183636.0_fp * phip2 +
                               5134574.0_fp * phip3 - 880548.0_fp * phip4) +
                      phi * (19365967.0_fp * phi - 65224244.0_fp * phip1 +
                             55053752.0_fp * phip2 - 23510468.0_fp * phip3 +
                             4067018.0_fp * phip4) +
                      phip1 * (56662212.0_fp * phip1 - 97838784.0_fp * phip2 +
                               42405032.0_fp * phip3 - 7408908.0_fp * phip4) +
                      phip2 * (43093692.0_fp * phip2 - 37913324.0_fp * phip3 +
                               6694608.0_fp * phip4) +
                      phip3 * (8449957.0_fp * phip3 - 3015728 * phip4) +
                      271779.0_fp * phip4 * phip4);

  const real beta5 = (phi * (6150211.0_fp * phi - 47460464.0_fp * phip1 +
                             76206736.0_fp * phip2 - 63394124.0_fp * phip3 +
                             27060170.0_fp * phip4 - 4712740.0_fp * phip5) +
                      phip1 * (94851237.0_fp * phip1 - 311771244.0_fp * phip2 +
                               262901672.0_fp * phip3 - 113206788.0_fp * phip4 +
                               19834350.0_fp * phip5) +
                      phip2 * (260445372.0_fp * phip2 - 444003904.0_fp * phip3 +
                               192596472.0_fp * phip4 - 33918804.0_fp * phip5) +
                      phip3 * (190757572.0_fp * phip3 - 166461044.0_fp * phip4 +
                               29442256.0_fp * phip5) +
                      phip4 * (36480687.0_fp * phip4 - 12950184.0_fp * phip5) +
                      1152561.0_fp * phip5 * phip5);

  const real alpha0 =
      (1.0_fp / 462.0_fp) / ((beta0 + 1e-10_fp) * (beta0 + 1e-10_fp));
  const real alpha1 =
      (5.0_fp / 77.0_fp) / ((beta1 + 1e-10_fp) * (beta1 + 1e-10_fp));
  const real alpha2 =
      (25.0_fp / 77.0_fp) / ((beta2 + 1e-10_fp) * (beta2 + 1e-10_fp));
  const real alpha3 =
      (100.0_fp / 231.0_fp) / ((beta3 + 1e-10_fp) * (beta3 + 1e-10_fp));
  const real alpha4 =
      (25.0_fp / 154.0_fp) / ((beta4 + 1e-10_fp) * (beta4 + 1e-10_fp));
  const real alpha5 =
      (1.0_fp / 77.0_fp) / ((beta5 + 1e-10_fp) * (beta5 + 1e-10_fp));

  const real alpha_sum_inv =
      1.0_fp / (alpha0 + alpha1 + alpha2 + alpha3 + alpha4 + alpha5);

  const real w0 = alpha0 * alpha_sum_inv;
  const real w1 = alpha1 * alpha_sum_inv;
  const real w2 = alpha2 * alpha_sum_inv;
  const real w3 = alpha3 * alpha_sum_inv;
  const real w4 = alpha4 * alpha_sum_inv;
  const real w5 = alpha5 * alpha_sum_inv;

  return w0 * p0 + w1 * p1 + w2 * p2 + w3 * p3 + w4 * p4 + w5 * p5;
};

template <uint ndofs, uint nd>
void YAKL_INLINE weno(SArray<real, 3, ndofs, nd, 2> &edgerecon,
                      SArray<real, 3, ndofs, nd, 1> const &dens) {
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      edgerecon(l, d, 0) = dens(l, d, 0);
      edgerecon(l, d, 1) = dens(l, d, 0);
    }
  }
}

template <uint ndofs, uint nd>
void YAKL_INLINE weno(SArray<real, 3, ndofs, nd, 2> &edgerecon,
                      SArray<real, 3, ndofs, nd, 3> const &dens) {
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      edgerecon(l, d, 0) =
          interp_weno3(dens(l, d, 2), dens(l, d, 1), dens(l, d, 0));
      edgerecon(l, d, 1) =
          interp_weno3(dens(l, d, 0), dens(l, d, 1), dens(l, d, 2));
    }
  }
}

template <uint ndofs, uint nd>
void YAKL_INLINE weno(SArray<real, 3, ndofs, nd, 2> &edgerecon,
                      SArray<real, 3, ndofs, nd, 5> const &dens) {
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      edgerecon(l, d, 0) =
          interp_weno5(dens(l, d, 4), dens(l, d, 3), dens(l, d, 2),
                       dens(l, d, 1), dens(l, d, 0));
      edgerecon(l, d, 1) =
          interp_weno5(dens(l, d, 0), dens(l, d, 1), dens(l, d, 2),
                       dens(l, d, 3), dens(l, d, 4));
    }
  }
}

template <uint ndofs, uint nd>
void YAKL_INLINE weno(SArray<real, 3, ndofs, nd, 2> &edgerecon,
                      SArray<real, 3, ndofs, nd, 7> const &dens) {
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      edgerecon(l, d, 0) = interp_weno7(
          dens(l, d, 6), dens(l, d, 5), dens(l, d, 4), dens(l, d, 3),
          dens(l, d, 2), dens(l, d, 1), dens(l, d, 0));
      edgerecon(l, d, 1) = interp_weno7(
          dens(l, d, 0), dens(l, d, 1), dens(l, d, 2), dens(l, d, 3),
          dens(l, d, 4), dens(l, d, 5), dens(l, d, 6));
    }
  }
}

template <uint ndofs, uint nd>
void YAKL_INLINE weno(SArray<real, 3, ndofs, nd, 2> &edgerecon,
                      SArray<real, 3, ndofs, nd, 9> const &dens) {
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      edgerecon(l, d, 0) =
          interp_weno9(dens(l, d, 8), dens(l, d, 7), dens(l, d, 6),
                       dens(l, d, 5), dens(l, d, 4), dens(l, d, 3),
                       dens(l, d, 2), dens(l, d, 1), dens(l, d, 0));
      edgerecon(l, d, 1) =
          interp_weno9(dens(l, d, 0), dens(l, d, 1), dens(l, d, 2),
                       dens(l, d, 3), dens(l, d, 4), dens(l, d, 5),
                       dens(l, d, 6), dens(l, d, 7), dens(l, d, 8));
    }
  }
}

template <uint ndofs, uint nd>
void YAKL_INLINE weno(SArray<real, 3, ndofs, nd, 2> &edgerecon,
                      SArray<real, 3, ndofs, nd, 11> const &dens) {
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      edgerecon(l, d, 0) = interp_weno11(
          dens(l, d, 10), dens(l, d, 9), dens(l, d, 8), dens(l, d, 7),
          dens(l, d, 6), dens(l, d, 5), dens(l, d, 4), dens(l, d, 3),
          dens(l, d, 2), dens(l, d, 1), dens(l, d, 0));
      edgerecon(l, d, 1) = interp_weno11(
          dens(l, d, 0), dens(l, d, 1), dens(l, d, 2), dens(l, d, 3),
          dens(l, d, 4), dens(l, d, 5), dens(l, d, 6), dens(l, d, 7),
          dens(l, d, 8), dens(l, d, 9), dens(l, d, 10));
    }
  }
}
} // namespace pamc
