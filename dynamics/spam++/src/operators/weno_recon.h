#ifndef _WENO_RECON_H_
#define _WENO_RECON_H_

#include "common.h"

real YAKL_INLINE interp_weno3(real phim1, real phi, real phip1){
    const real p0 = (-1.0/2.0) * phim1 + (3.0/2.0) * phi;
    const real p1 = (1.0/2.0) * phi + (1.0/2.0) * phip1;

    const real beta1 = (phip1 - phi) * (phip1 - phi);
    const real beta0 = (phi - phim1) * (phi - phim1);

    const real alpha0 = (1.0/3.0) /((beta0 + 1e-10) * (beta0 + 1.0e-10));
    const real alpha1 = (2.0/3.0)/((beta1 + 1e-10) * (beta1 + 1.0e-10));

    const real alpha_sum_inv = 1.0/(alpha0 + alpha1);

    const real w0 = alpha0 * alpha_sum_inv;
    const real w1 = alpha1 * alpha_sum_inv;

    return w0 * p0 + w1 * p1;
};

real YAKL_INLINE interp_weno5(real phim2, real phim1, real phi, real phip1, real phip2){

   const real p0 = (1.0/3.0)*phim2 - (7.0/6.0)*phim1 + (11.0/6.0)*phi;
   const real p1 = (-1.0/6.0) * phim1 + (5.0/6.0)*phi + (1.0/3.0)*phip1;
   const real p2 = (1.0/3.0) * phi + (5.0/6.0)*phip1 - (1.0/6.0)*phip2;

   const real beta2 = (13.0/12.0 * (phi - 2.0 * phip1 + phip2)*(phi - 2.0 * phip1 + phip2)
                        + 0.25 * (3.0 * phi - 4.0 * phip1 + phip2)*(3.0 * phi - 4.0 * phip1 + phip2));
   const real beta1 = (13.0/12.0 * (phim1 - 2.0 * phi + phip1)*(phim1 - 2.0 * phi + phip1)
                        + 0.25 * (phim1 - phip1)*(phim1 - phip1));
   const real beta0 = (13.0/12.0 * (phim2 - 2.0 * phim1 + phi)*(phim2 - 2.0 * phim1 + phi)
                        + 0.25 * (phim2 - 4.0 * phim1 + 3.0 * phi)*(phim2 - 4.0 * phim1 + 3.0 * phi));

   const real alpha0 = 0.1/((beta0 + 1e-10) * (beta0 + 1e-10));
   const real alpha1 = 0.6/((beta1 + 1e-10) * (beta1 + 1e-10));
   const real alpha2 = 0.3/((beta2 + 1e-10) * (beta2 + 1e-10));

   const real alpha_sum_inv = 1.0/(alpha0 + alpha1 + alpha2);

   const real w0 = alpha0 * alpha_sum_inv;
   const real w1 = alpha1 * alpha_sum_inv;
   const real w2 = alpha2 * alpha_sum_inv;


   return w0 * p0 + w1 * p1 + w2 * p2;
};

real YAKL_INLINE interp_weno7(real phim3, real phim2, real phim1, real phi, real phip1, real phip2, real phip3){

    const real p0 = (-1.0/4.0)*phim3 + (13.0/12.0) * phim2 + (-23.0/12.0) * phim1 + (25.0/12.0)*phi;
    const real p1 = (1.0/12.0)*phim2 + (-5.0/12.0)*phim1 + (13.0/12.0)*phi + (1.0/4.0)*phip1;
    const real p2 = (-1.0/12.0)*phim1 + (7.0/12.0)*phi +  (7.0/12.0)*phip1 + (-1.0/12.0)*phip2;
    const real p3 = (1.0/4.0)*phi + (13.0/12.0)*phip1 + (-5.0/12.0)*phip2 + (1.0/12.0)*phip3;


    const real beta0 = (phim3*(547.0*phim3 - 3882.0*phim2 + 4642.0*phim1 - 1854.0*phi)
                         + phim2*(7043.0*phim2 - 17246.0*phim1 + 7042.0*phi)
                         + phim1*(11003.0*phim1 - 9402.0*phi)
                         + 2107.0*phi*phi);
    const real beta1 =(phim2*(267.0*phim2 - 1642.0*phim1 + 1602.0*phi - 494.0*phip1)
                        + phim1*(2843.0*phim1 - 5966.0*phi + 1922.0*phip1)
                        + phi*(3443.0*phi - 2522.0*phip1)
                        + 547.0*phip1*phip1);
    const real beta2 = (phim1*(547.0*phim1 - 2522.0*phi + 1922.0*phip1 - 494.0*phip2)
                         + phi*(3443.0*phi -5966.0*phip1 + 1602.0*phip2)
                         + phip1*(2843.0*phip1 - 1642.0*phip2)
                         + 267.0*phip2* phip2);
    const real beta3 = (phi*(2107.0*phi - 9402.0*phip1 + 7042.0*phip2 - 1854.0*phip3)
                         + phip1*(11003.0*phip1 - 17246.0*phip2 + 4642.0*phip3)
                         + phip2*(7043.0*phip2 - 3882.0*phip3)
                         + 547.0*phip3*phip3);

    const real alpha0 = (1.0/35.0)/((beta0 + 1e-10) * (beta0 + 1e-10));
    const real alpha1 = (12.0/35.0)/((beta1 + 1e-10) * (beta1 + 1e-10));
    const real alpha2 = (18.0/35.0)/((beta2 + 1e-10) * (beta2 + 1e-10));
    const real alpha3 = (4.0/35.0)/((beta3 + 1e-10) * (beta3 + 1e-10));

    const real alpha_sum_inv = 1.0/(alpha0 + alpha1 + alpha2 + alpha3);

    const real w0 = alpha0 * alpha_sum_inv;
    const real w1 = alpha1 * alpha_sum_inv;
    const real w2 = alpha2 * alpha_sum_inv;
    const real w3 = alpha3 * alpha_sum_inv;


    return w0 * p0 + w1 * p1 + w2 * p2 + w3 * p3;
};

real YAKL_INLINE interp_weno9(real phim4, real phim3, real phim2, real phim1, real phi, real phip1, real phip2, real phip3, real phip4){
    const real p0 = (1.0/5.0) * phim4 + (-21.0/20.0)*phim3 +  (137.0/60.0)*phim2 +  (-163.0/60.0)*phim1 +  (137.0/60.0)*phi;
    const real p1 = (-1.0/20.0)*phim3 + (17.0/60.0)*phim2 + (-43.0/60.0)*phim1 + (77.0/60.0)*phi + (1.0/5.0)*phip1;
    const real p2 = (1.0/30.0)*phim2 + (-13.0/60.0)*phim1 + (47.0/60.0)*phi +  (9.0/20.0)*phip1 + (-1.0/20.0)*phip2;
    const real p3 = (-1.0/20.0)*phim1 + (9.0/20.0)*phi + (47.0/60.0)*phip1 + (-13.0/60.0)*phip2 + (1.0/30.0)*phip3;
    const real p4 = (1.0/5.0)*phi + (77.0/60.0)*phip1 + (-43.0/60.0)*phip2 + (17.0/60.0)*phip3 + (-1.0/20.0)*phip4;

    const real beta0 = (phim4 * (22658.0 * phim4 - 208501.0 * phim3 + 364863.0 * phim2 - 288007.0 * phim1 + 86329.0 * phi)
                         + phim3 * (482963.0 * phim3  - 1704396.0 * phim2 + 1358458.0 * phim1 - 411487.0 * phi)
                         + phim2 * (1521393.0 * phim2 - 2462076.0 * phim1 + 758823.0 * phi )
                         + phim1 * (1020563.0 * phim1 - 649501.0 * phi)
                         + 107918.0 * phi * phi );
    const real beta1 = (phim3 * (6908.0*phim3 - 60871.0 * phim2 + 99213.0*phim1 - 70237.0 * phi + 18079.0 * phip1)
                         +phim2 * (138563.0 * phim2 - 464976.0 * phim1 + 337018.0*phi - 88297.0 * phip1)
                         + phim1 * (406293.0 * phim1 - 611976.0 * phi + 165153.0 * phip1)
                         + phi * (242723.0 * phi - 140251.0 * phip1)
                         +22658.0 * phip1 * phip1);
    const real beta2 = (phim2 * (6908.0 * phim2 - 51001.0 * phim1 + 67923.0 * phi - 38947.0 * phip1 + 8209.0 * phip2)
                         +phim1 * (104963.0 * phim1 - 299076.0 * phi + 179098.0 * phip1 - 38947.0 * phip2)
                         +phi * (231153.0 * phi - 299076.0 * phip1 + 67923.0*phip2)
                         +phip1 * (104963.0 * phip1 - 51001.0 * phip2)
                         + 6908.0 * phip2 * phip2);
    const real beta3 = (phim1 * (22658.0 * phim1 - 140251.0 * phi + 165153.0 * phip1 - 88297.0 * phip2 + 18079.0 * phip3)
                         + phi * (242723.0 * phi - 611976.0 * phip1 + 337018.0 * phip2 - 70237.0 * phip3)
                         +phip1 * (406293.0 * phip1 - 464976.0 * phip2 + 99213.0 * phip3 )
                         +phip2 * (138563.0 * phip2 - 60871.0 * phip3)
                         + 6908.0 * phip3 * phip3);
    const real beta4 = (phi * (107918.0 * phi - 649501.0 * phip1 + 758823.0 * phip2 - 411487.0 * phip3 + 86329.0 * phip4)
                         +phip1 * (1020563.0 * phip1 - 2462076.0 * phip2 + 1358458.0 * phip3 - 288007.0 * phip4)
                         +phip2 * (1521393.0 * phip2 - 1704396.0 * phip3 + 364863.0*phip4)
                         +phip3 * (482963.0 * phip3 - 208501.0 * phip4)
                         + 22658.0 * phip4 * phip4 );


    const real alpha0 = (1.0/126.0)/((beta0 + 1e-10)*(beta0 + 1e-10));
    const real alpha1 = (10.0/63.0)/((beta1 + 1e-10)*(beta1 + 1e-10));
    const real alpha2 = (10.0/21.0)/((beta2 + 1e-10)*(beta2 + 1e-10) );
    const real alpha3 = (20.0/63.0)/((beta3 + 1e-10)*(beta3 + 1e-10));
    const real alpha4 = (5.0/126.0)/((beta4 + 1e-10)*(beta4 + 1e-10));


    const real alpha_sum_inv = 1.0/(alpha0 + alpha1 + alpha2 + alpha3 + alpha4 );


    const real w0 = alpha0*alpha_sum_inv;
    const real w1 = alpha1*alpha_sum_inv;
    const real w2 = alpha2*alpha_sum_inv;
    const real w3 = alpha3*alpha_sum_inv;
    const real w4 = alpha4*alpha_sum_inv;

    return w0 * p0 + w1 * p1 + w2 * p2 + w3 * p3 + w4 * p4;
};

real YAKL_INLINE interp_weno11(real phim5, real phim4, real phim3, real phim2, real phim1, real phi, real phip1, real phip2, real phip3, real phip4, real phip5  ){

    const real p0 = ((-1.0/6.0) * phim5 + (31.0/30.0)*phim4 +  (-163.0/60.0)*phim3 +  (79.0/20.0)*phim2
                     +  (-71.0/20.0)*phim1 +  (49.0/20.0)*phi);
    const real p1 = ((1.0/30.0)*phim4 + (-13.0/60.0)*phim3 + (37.0/60.0)*phim2 + (-21.0/20.0)*phim1
                     + (29.0/20.0)*phi + (1.0/6.0)*phip1);
    const real p2 = ((-1.0/60.0)*phim3 + (7.0/60.0)*phim2 + (-23.0/60.0)*phim1 +  (19.0/20.0)*phi
                     + (11.0/30.0)*phip1 + (-1.0/30.0)*phip2);
    const real p3 = ((1.0/60.0)*phim2 + (-2.0/15.0)*phim1 + (37.0/60.0)*phi + (37.0/60.0)*phip1
                   + (-2.0/15.0)*phip2 + (1.0/60.0)*phip3);
    const real p4 = ((-1.0/30.0)*phim1 + (11.0/30.0)*phi + (19.0/20.0)*phip1 + (-23.0/60.0)*phip2
                     + (7.0/60.0)*phip3 + (-1.0/60.0)*phip4);
    const real p5 = ((1.0/6.0)*phi + (29.0/20.0)*phip1 + (-21.0/20.0)*phip2 + (37.0/60.0)*phip3
                     + (-13.0/60.0)*phip4 + (1.0/30.0)*phip5);


    const real beta0 = ( phim5 * (1152561.0*phim5 - 12950184.0*phim4 + 29442256.0*phim3 - 33918804.0*phim2 + 19834350.0*phim1 -4712740.0*phi)
                         + phim4 * (36480687.0*phim4 - 166461044.0*phim3 + 192596472.0*phim2 - 113206788.0*phim1 + 27060170.0*phi)
                         + phim3 * (190757572.0*phim3 -  444003904.0*phim2 + 262901672.0*phim1 - 63394124.0*phi)
                         +phim2*(260445372.0*phim2 - 311771244.0*phim1 + 76206736.0*phi)
                         +phim1*(94851237.0*phim1-47460464.0*phi)
                         +6150211.0*phi*phi);

    const real beta1 = (phim4 * (271779.0*phim4 - 3015728.0*phim3 + 6694608.0*phim2 - 7408908.0*phim1 + 4067018.0*phi - 880548.0*phip1)
                         + phim3 *(8449957.0*phim3 - 37913324.0*phim2 + 42405032.0*phim1 - 23510468.0*phi + 5134574.0*phip1 )
                         + phim2 *(43093692.0*phim2 -  97838784.0*phim1 + 55053752.0*phi - 12183636.0*phip1)
                         + phim1 *(56662212.0*phim1 - 65224244.0*phi + 14742480.0*phip1)
                         + phi * (19365967.0*phi-9117992.0*phip1)
                         + 1152561.0*phip1*phip1);

    const real beta2 =  (phim3 * (139633.0*phim3 - 1429976.0*phim2 + 2863984.0*phim1 - 2792660.0*phi + 1325006.0*phip1 - 245620.0*phip2)
                          + phim2*(3824847.0*phim2 - 15880404.0*phim1 + 15929912.0*phi - 7727988.0*phip1 + 1458762.0*phip2)
                          + phim1*(17195652.0*phim1 - 35817664.0*phi + 17905032.0*phip1 - 3462252.0*phip2)
                          + phi*(19510972.0*phi - 20427884.0*phip1 + 4086352.0*phip2)
                          + phip1*(5653317.0*phip1 - 2380800.0*phip2)
                          + 271779.0*phip2*phip2);

    const real beta3 = (phim2*(271779.0*phim2 - 2380800.0*phim1 + 4086352.0*phi - 3462252.0*phip1 + 1458762.0*phip2 - 245620.0*phip3)
                         + phim1*(5653317.0*phim1 - 20427884.0*phi + 17905032.0*phip1 - 7727988.0*phip2 + 1325006.0*phip3)
                         + phi*(19510972.0*phi - 35817664.0*phip1 + 15929912.0*phip2 - 2792660.0*phip3)
                         + phip1*(17195652.0*phip1 - 15880404.0*phip2 + 2863984.0*phip3)
                         + phip2*(3824847.0*phip2 - 1429976.0*phip3)
                         + 139633.0*phip3*phip3);

    const real beta4 = (phim1*(1152561.0*phim1 - 9117992.0*phi + 14742480.0*phip1 - 12183636.0*phip2 + 5134574.0*phip3 - 880548.0*phip4)
                         + phi*(19365967.0*phi - 65224244.0*phip1 + 55053752.0*phip2 - 23510468.0*phip3 + 4067018.0*phip4)
                         + phip1*(56662212.0*phip1 - 97838784.0*phip2 + 42405032.0*phip3 - 7408908.0*phip4)
                         + phip2*(43093692.0*phip2 - 37913324.0*phip3 + 6694608.0*phip4)
                         + phip3*(8449957.0*phip3 - 3015728*phip4)
                         + 271779.0*phip4*phip4);

    const real beta5 = (phi*(6150211.0*phi - 47460464.0*phip1 + 76206736.0*phip2 - 63394124.0*phip3 + 27060170.0*phip4 - 4712740.0*phip5)
                         + phip1*(94851237.0*phip1 - 311771244.0*phip2 + 262901672.0*phip3 - 113206788.0*phip4 + 19834350.0*phip5)
                         + phip2*(260445372.0*phip2 - 444003904.0*phip3 + 192596472.0*phip4 - 33918804.0*phip5)
                         + phip3*(190757572.0*phip3 - 166461044.0*phip4 + 29442256.0*phip5)
                         + phip4*(36480687.0*phip4 - 12950184.0*phip5)
                         + 1152561.0*phip5 * phip5);


    const real alpha0 = (1.0/462.0)/((beta0 + 1e-10) * (beta0 + 1e-10));
    const real alpha1 = (5.0/77.0)/((beta1 + 1e-10) * (beta1 + 1e-10));
    const real alpha2 = (25.0/77.0)/((beta2 + 1e-10) * (beta2 + 1e-10));
    const real alpha3 = (100.0/231.0)/((beta3 + 1e-10) * (beta3 + 1e-10));
    const real alpha4 = (25.0/154.0)/((beta4 + 1e-10) * (beta4 + 1e-10));
    const real alpha5 = (1.0/77.0)/((beta5 + 1e-10) * (beta5 + 1e-10));

    const real alpha_sum_inv = 1.0/(alpha0 + alpha1 + alpha2 + alpha3 + alpha4 + alpha5);

    const real w0 = alpha0*alpha_sum_inv;
    const real w1 = alpha1*alpha_sum_inv;
    const real w2 = alpha2*alpha_sum_inv;
    const real w3 = alpha3*alpha_sum_inv;
    const real w4 = alpha4*alpha_sum_inv;
    const real w5 = alpha5*alpha_sum_inv;

    return w0 * p0 + w1 * p1 + w2 * p2 + w3 * p3 + w4 * p4 + w5 * p5;
};





template<uint ndofs, uint nd> void YAKL_INLINE weno(SArray<real,ndofs,nd,2> &edgerecon, SArray<real,ndofs,nd,1> const &dens) {
    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {
          edgerecon(l,d,0) = dens(l,d,0);
          edgerecon(l,d,1) = dens(l,d,0);
        }
      }
}

template<uint ndofs, uint nd> void YAKL_INLINE weno(SArray<real,ndofs,nd,2> &edgerecon, SArray<real,ndofs,nd,3> const &dens) {
    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {
          edgerecon(l,d,0) = interp_weno3(dens(l,d,2), dens(l,d,1), dens(l,d,0));
          edgerecon(l,d,1) = interp_weno3(dens(l,d,0), dens(l,d,1), dens(l,d,2));
        }
      }
}

template<uint ndofs, uint nd> void YAKL_INLINE weno(SArray<real,ndofs,nd,2> &edgerecon, SArray<real,ndofs,nd,5> const &dens) {
    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {
          edgerecon(l,d,0) = interp_weno5(dens(l,d,4), dens(l,d,3), dens(l,d,2), dens(l,d,1), dens(l,d,0));
          edgerecon(l,d,1) = interp_weno5(dens(l,d,0), dens(l,d,1), dens(l,d,2), dens(l,d,3), dens(l,d,4));
        }
      }
}

template<uint ndofs, uint nd> void YAKL_INLINE weno(SArray<real,ndofs,nd,2> &edgerecon, SArray<real,ndofs,nd,7> const &dens) {
    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {
          edgerecon(l,d,0) = interp_weno7(dens(l,d,6), dens(l,d,5), dens(l,d,4), dens(l,d,3), dens(l,d,2), dens(l,d,1), dens(l,d,0));
          edgerecon(l,d,1) = interp_weno7(dens(l,d,0), dens(l,d,1), dens(l,d,2), dens(l,d,3), dens(l,d,4), dens(l,d,5), dens(l,d,6));
        }
      }
}

template<uint ndofs, uint nd> void YAKL_INLINE weno(SArray<real,ndofs,nd,2> &edgerecon, SArray<real,ndofs,nd,9> const &dens) {
    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {
          edgerecon(l,d,0) = interp_weno9(dens(l,d,8), dens(l,d,7), dens(l,d,6), dens(l,d,5), dens(l,d,4), dens(l,d,3), dens(l,d,2), dens(l,d,1), dens(l,d,0));
          edgerecon(l,d,1) = interp_weno9(dens(l,d,0), dens(l,d,1), dens(l,d,2), dens(l,d,3), dens(l,d,4), dens(l,d,5), dens(l,d,6), dens(l,d,7), dens(l,d,8));
        }
      }
}

template<uint ndofs, uint nd> void YAKL_INLINE weno(SArray<real,ndofs,nd,2> &edgerecon, SArray<real,ndofs,nd,11> const &dens) {
    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {
          edgerecon(l,d,0) = interp_weno11(dens(l,d,10), dens(l,d,9), dens(l,d,8), dens(l,d,7), dens(l,d,6), dens(l,d,5), dens(l,d,4), dens(l,d,3), dens(l,d,2), dens(l,d,1), dens(l,d,0));
          edgerecon(l,d,1) = interp_weno11(dens(l,d,0), dens(l,d,1), dens(l,d,2), dens(l,d,3), dens(l,d,4), dens(l,d,5), dens(l,d,6), dens(l,d,7), dens(l,d,8), dens(l,d,9), dens(l,d,10));
        }
      }
}
#endif
