
#include "const.h"


/*
  Gives a cosine ellipsiod centered at (x0,y0,z0) with radius (xrad,yrad,zrad) and amplitude amp
*/
YAKL_INLINE real ellipsoid_cosine(real x   , real y   , real z ,
                                  real x0  , real y0  , real z0,
                                  real xrad, real yrad, real zrad, real amp, real pwr) {
  real val = 0;
  real xn = (x-x0)/xrad;
  real yn = (y-y0)/yrad;
  real zn = (z-z0)/zrad;
  real dist = sqrt( xn*xn + yn*yn + zn*zn );
  if (dist <= 1._fp) {
    val = amp * pow( (cos(M_PI*dist)+1)/2 , pwr );
  }
  return val;
}



int main() {
  int constexpr hs = 1;
  int constexpr nz = 200;
  int constexpr nx = 200;
  real constexpr dx = 0.2;
  realHost2d rho_u("rho_u",nz+2*hs,nx+2*hs);
  realHost2d rho_w("rho_w",nz+2*hs,nx+2*hs);

  // Load random values
  for (int k=0; k < nz; k++) {
    for (int i=0; i < nx; i++) {
      rho_u(hs+k,hs+i) = ellipsoid_cosine(i,1,k,nx/2,1,nz/2,nx/4,1,nz/4,2,2) + 2;
      rho_w(hs+k,hs+i) = ellipsoid_cosine(i,1,k,nx/2,1,nz/2,nx/4,1,nz/4,2,10);
    }
  }

  // Set z-direction ghost cells
  for (int i=0; i < nx+2*hs; i++) {
    rho_u(0,i) = rho_u(1,i);
    rho_w(0,i) = 0;
    rho_u(nz+1,i) = rho_u(nz,i);
    rho_w(nz+1,i) = 0;
  }

  // Set x-direction ghost cells
  for (int k=0; k < nz+2*hs; k++) {
    rho_u(k,0) = rho_u(k,nx);
    rho_w(k,0) = rho_w(k,nx);
    rho_u(k,nx+1) = rho_u(k,1);
    rho_w(k,nx+1) = rho_w(k,1);
  }

  realHost2d phix("phix",nz,nx);
  realHost2d phiz("phiz",nz,nx);

  memset(phix,0._fp);
  memset(phiz,0._fp);

  for (int k=0; k < nz; k++) {
    phix(k,0) = (rho_u(hs+k,hs+0+1) - rho_u(hs+k,hs+0-1))/2;
    real avg = phix(k,0);
    for (int i=1; i < nx; i++) {
      phix(k,i) = phix(k,i-1) + (rho_u(hs+k,hs+i+1) - rho_u(hs+k,hs+i-1))/2;
      avg += phix(k,i);
    }
    avg /= nx;
    for (int i=0; i < nx; i++) {
      phix(k,i) -= avg;
    }
  }

  for (int i=0; i < nx; i++) {
    phiz(0,i) = (rho_w(hs+0+1,hs+i) - rho_w(hs+0-1,hs+i))/2;
    real avg = phiz(0,i);
    for (int k=1; k < nz; k++) {
      phiz(k,i) = phiz(k-1,i) + (rho_w(hs+k+1,hs+i) - rho_w(hs+k-1,hs+i))/2;
      avg += phiz(k,i);
    }
    avg /= nz;
    for (int k=0; k < nz; k++) {
      phiz(k,i) -= avg;
    }
  }


  for (int k=0; k < nz; k++) {
    phix(k,0) = (rho_u(hs+k,hs+0+1) - rho_u(hs+k,hs+0-1))/2;
    real avg = phix(k,0);
    for (int i=1; i < nx; i++) {
      phix(k,i) = phix(k,i-1) + (rho_u(hs+k,hs+i+1) - rho_u(hs+k,hs+i-1))/2;
      avg += phix(k,i);
    }
    avg /= nx;
    for (int i=0; i < nx; i++) {
      phix(k,i) -= avg;
    }
  }

  realHost2d phi("phi",nz+2*hs,nx+2*hs);
  memset(phi,0._fp);

  real sum = 0;
  for (int k=0; k < nz; k++) {
    for (int i=0; i < nx; i++) {
      phi(hs+k,hs+i) = phix(k,i) + phiz(k,i);
      // phi(hs+k,hs+i) = phix(k,i);
      // phi(hs+k,hs+i) = phiz(k,i);
      sum += phi(hs+k,hs+i);
    }
  }

  std::cout << sum << "\n";

  // // Set z-direction ghost cells
  // for (int i=0; i < nx+2*hs; i++) {
  //   phi(0,i) = 0;
  //   phi(nz+1,i) = 0;
  // }

  // // Set x-direction ghost cells
  // for (int k=0; k < nz+2*hs; k++) {
  //   phi(k,0) = phi(k,nx);
  //   phi(k,nx+1) = phi(k,1);
  // }


  real norm = 0;
  real norm_denom = 0;
  for (int k=1; k < nz-1; k++) {
    for (int i=1; i < nx-1; i++) {
      real phi_grad = ( phi  (hs+k+1,hs+i) - phi  (hs+k-1,hs+i) ) / 2 + 
                      ( phi  (hs+k,hs+i+1) - phi  (hs+k,hs+i-1) ) / 2 ;
      real mom_grad = ( rho_w(hs+k+1,hs+i) - rho_w(hs+k-1,hs+i) ) / 2 + 
                      ( rho_u(hs+k,hs+i+1) - rho_u(hs+k,hs+i-1) ) / 2 ;

      // real phi_grad = ( phi  (hs+k,hs+i+1) - phi  (hs+k,hs+i-1) ) / 2 ;
      // real mom_grad = ( rho_u(hs+k,hs+i+1) - rho_u(hs+k,hs+i-1) ) / 2 ;

      // real phi_grad = ( phi  (hs+k+1,hs+i) - phi  (hs+k-1,hs+i) ) / 2 ;
      // real mom_grad = ( rho_w(hs+k+1,hs+i) - rho_w(hs+k-1,hs+i) ) / 2 ; 

      norm += abs(mom_grad - phi_grad);
      norm_denom += abs(mom_grad);
    }
  }

  realHost2d mag("mag",nz,nx);

  std::cout << norm / norm_denom << "\n";
  
}

