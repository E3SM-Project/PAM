#ifndef _HAMILTONIAN_H_
#define _HAMILTONIAN_H_

#include "common.h"

// this defines Hamiltonian, and associated functional derivatives
// It is H = K + I + P
// K depends on slice vs. no-slice, and split vs. unsplit
// I + P are function of thermodynamics + entropic variable we use, and also equation set


// ie we have H[v,w,vt,density,densityfct,p] + drop w,vt,p as needed
//K is function of total_density, v,w,vt
//I,P are functions are density,densityfct,p

//actually I can write something SUPER general here, with some conversion functions that do nothing for most cases...
// really difference is in velocity choice and p/no-p

// for AN/FC, go from density + densityfct vars to total density (=alpha), entropic variable (=entropic variable density/total_density), concentrations

// basically we can predict rho or rho_d 
// p vs. no p?
// can write Hamiltonian + FD in a general way for ANY internal energy/entropic variable, for rho/rho_d and p/no-p variants

// have functions like dH/d(density), dH/d(densityFCT), etc.
// this fits the GENERAL hamiltonian formulation we are developing

// also need code to get total density from density + densityfct vars...
// this is not too difficult!

// then convert predicted variables to alpha, entropic variable, concentrations 

// can probably merge swe/tswe stuff into here as well...



// 4 types of Hk: v, v+w, v+w+vt, v+vt
// 2 types of Hs (=Hp + Hi, sort of like static energy): D, D+p : these are also specialized to a given eqn set and variables ie rho vs. rho_d
//ie split/unsplit, slice/no-slice, p/no-p

// Maybe specific Hamiltonians inherit from 1 of each?
// Yes this is the way...


//Need a class or function that handles
// 1) getting h/D from dens + densfct
// 2) getting alpha, theta, qn + associated partial derivatives from dens + densfct
// 3) Setting initial conditions ie values for dens/densfct given h,s (swe/tswe) or p,T,qs (ce)
// This should then enable use of arbitrary linear combinations of component densities for ce; and also handle the swe/tswe case

class progset {
public:
  bool is_initialized;
  
  progset() {
    this->is_initialized = false;
  }

  void initialize()
  {
    this->is_initialized = true;
  }
  
  virtual real YAKL_INLINE compute_D(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k) {};
  //template<uint kdof> virtual real YAKL_INLINE compute_dDdDk(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k) {};
  virtual real YAKL_INLINE compute_alpha(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k) {};
  virtual real YAKL_INLINE compute_theta(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k) {};
  //template<uint n> virtual real YAKL_INLINE compute_qn(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k) {};
  // INITIAL CONDITION SETTING?
  // VARIOUS PARTIAL DERIVATIVES?
};  

class hvS : public progset {
  real YAKL_INLINE compute_D(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
    {return dens(0, k+ks, j+js, i+is);}
};

class rhovS : public progset {   //single component
  real YAKL_INLINE compute_D(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
    {return dens(0, k+ks, j+js, i+is);}
}; 
class rhonvS : public progset { //multicomponent rho, rho_s
  real YAKL_INLINE compute_D(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
    {return dens(0, k+ks, j+js, i+is);}
}; 
class rhodvS : public progset { //multicomponent rho_d, rho_s (s=v,l,i)
  real YAKL_INLINE compute_D(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
    {return dens(0, k+ks, j+js, i+is) + densfct(0, k+ks, j+js, i+is) + densfct(1, k+ks, j+js, i+is) + densfct(2, k+ks, j+js, i+is);}
}; 




class pvpe_2D {
public:
  const progset *xv;
  bool is_initialized;

   pvpe_2D() {
     this->is_initialized = false;
}

void initialize(const progset &xvar)
{
  this->xv = &xvar;
  this->is_initialized = true;
}

void YAKL_INLINE compute_q0f0(realArr q0, realArr f0, const realArr v, const realArr dens, const realArr densfct, const realArr coriolis, int is, int js, int ks, int i, int j, int k)
{
  real hv;
  SArray<real,4> D;

  // compute zeta = D2 v
  //compute_D2<1>(q0, v, is, js, ks, i, j, k);

  // compute q0 = zeta / R h
  //D(0) = xvar.compute_D(dens, densfct, is, js, ks, i, j, k);
  //ADD/FIX THE REST OF THIS...
    //compute_R<0> (hv, D, is, js, ks, i, j, k);
    
    //q0(0, k+ks, j+js, i+is) = q0(0, k+ks, j+js, i+is) / hv;
    //f0(0, k+ks, j+js, i+is) = coriolis(0, k+ks, j+js, i+is) / hv;

}

void YAKL_INLINE compute_PVPE(realArr v, realArr dens, realArr densfct, int is, int js, int ks, int i, int j, int k)
{
//K = 1/2 * h * v * v;  
}

};








class v_Hk {
  
public:
  const progset *xv;
  Geometry<ndims,1,1,1> *primal_geometry;
  Geometry<ndims,1,1,1> *dual_geometry;
  bool is_initialized;

   v_Hk() {
     this->is_initialized = false;
}

void initialize(const progset &xvar, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
{
  this->xv = &xvar;
  this->primal_geometry = &primal_geom;
  this->dual_geometry = &dual_geom;
  this->is_initialized = true;
}

void YAKL_INLINE compute_KE(const realArr v, const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
{
  //h = xvar.compute_D(densvar, densfctvar, is, js, ks, i, j, k);
//K = 1/2 * h * v * v;  
}


 void YAKL_INLINE compute_dKdx(realArr F, realArr K, realArr HE, const realArr v, const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k)
{
  //Sarray<real,2,2> U;
  //Sarray<real,2,2> vcell;
  //Sarray<real,2> D0;

  //Compute U = H * v
//THIS ONE DOESN'T EXIST/IS BROKEN...
  //compute_H_cell<1, diff_ord>(U, vcell, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

  // compute he = phi * h0
  // Exploits linearity of Hodge star to compute D0 = I D = I "\sum dens" = "\sum" I dens = "\sum" dens0
  //FIX THE THE D0(1) COMPUTATION
  //D0(0) = xvar.compute_D(dens0, densfct0, is, js, ks, i, j, k);
  //D0(1) = xvar.compute_D(dens0, densfct0, is, js, ks, i, j, k);
//THIS ONE DOESN'T EXIST/IS BROKEN...
  //compute_phi<0>(HE, D0, is, js, ks, i, j, k);

  //compute F = he * U
  // FIX INDEXING HERE...
  //Fvar(0, k+ks, j+js, i+is) = U(0) * HEvar(0, k+dks, j+js, i+is);
  //Fvar(1, k+ks, j+js, i+is) = U(1) * HEvar(1, k+dks, j+js, i+is);

  //compute_phiT(K, U, vcell, is, js, ks, i, j, k);
  //K(0, k+ks, j+js, i+is) *= 0.5;
  

}

 };
 
 
 // class vw_Hk {
 // public:
 //   bool is_initialized;
 // 
 //    vw_Hk() {
 //      this->is_initialized = false;
 // }
 // 
 // void YAKL_INLINE compute_H(const realArr v, const realArr w, const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 // {
 //   HK = 1/2 * h * v * v + 1/2 * h * w * w;
 // };
 // 
 //  void YAKL_INLINE compute_dHKdx(realArr F, realArr Fw, const realArr v, const realArr vw, const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 // {
 //   F = h * v;
 //   Fw = h * w;
 //   //CAREFUL WITH THIS! Specialize to rho vs. rho_d?
 //   // ADD OR REPLACE MODE?
 //   B = 1/2 * v * v + 1/2 w * w;
 //   Bfct = 1/2 * v * v;
 // };
 // 
 //  };

// template<uint nt, uint ntfct> class swe_Hs {
//   public:
//     bool is_initialized;
// 
//      swe_Hs() {
//        this->is_initialized = false;
//   }
// 
//   void YAKL_INLINE compute_P(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
//   {
//     P = g * h * hs + 1/2 g * h * h + sum_nt 1/2 h * t + sum_nt 1/2 h * tfct
//   };
// 
//   void YAKL_INLINE compute_I(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
//   {
//     I = 0
//   };
// 
//    void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
//    {
//      // ADD OR REPLACE MODE?
//      dHdh = g hs + gh + sum_nt 1/2 t + sum_nt 1/2 tfct;     
//      dHdt = 1/2 h;
//      dHdtfct = 1/2 h;
//    };
// 
//  };
  
  
template<uint nt, uint ntfct> class tswe_Hs {
 public:
   const progset *xv;
   Geometry<ndims,1,1,1> *primal_geometry;
   Geometry<ndims,1,1,1> *dual_geometry;
   bool is_initialized;
   realArr hs0;
   
    tswe_Hs() {
      this->is_initialized = false;
 }
 
 void initialize(const progset &xvar, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom) 
   //const realArr &hs)
 {
   this->xv = &xvar;
   this->primal_geometry = &primal_geom;
   this->dual_geometry = &dual_geom;
   //compute_I<1, diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   this->is_initialized = true;
 }
 
 void YAKL_INLINE compute_PE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 {
   //P = S * hs + 1/2 S * h + sum_nt 1/2 h * t + sum_nt 1/2 h * tfct;
 }

 void YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 {
   //I = 0;
 }

 //HOW SHOULD HS + OTHER CONSTANTS (GEOPOTENTIAL, ETC.) BE HANDLED IN A UNIFORM WAY
 //IDEALLY THEY ARE PART OF INITIALIZE FOR HAMILTONIAN 
 // MAYBE INITIALIZE TAKES CONSTANT VARS AS ARGUMENTS?
 // THEN WE CAN ALSO DO THINGS LIKE TAKE I HS ONCE, I GEOP ONCE, ETC.
// YES DO IT LIKE THIS...

  void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr hs, int is, int js, int ks, int i, int j, int k)
  {
    // dens0 is h,S,t; densfct is tfct


//ELIMINATE THIS EVENTUALLY...
    //compute I hs
    SArray<real,1> hs0;
    compute_I<1, diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    
    // Compute dHdh = 1/2 S + sum_nt 1/2 t + sum_nt 1/2 tfct
    B(0, k+ks, j+js, i+is) = dens0(1, k+ks, j+js, i+is)/2. ;
    for (int l=2; l<2+nt; l++) { B(0, k+ks, j+js, i+is) += dens0(l, k+ks, j+js, i+is)/2.;}
    for (int l=0; l<ntfct; l++) { B(0, k+ks, j+js, i+is) += densfct0(l, k+ks, j+js, i+is)/2.;}
          
    //Compute dHdS = 1/2 h + hs
    B(1, k+ks, j+js, i+is) = hs0(0) + dens0(0, k+ks, j+js, i+is)/2.;

    //Compute dHdt = 1/2 h
    for (int l=2; l<2+nt; l++) { B(l, k+ks, j+js, i+is) = dens0(0, k+ks, j+js, i+is)/2.;}
    
    //Compute dHdtfct = 1/2 h
    for (int l=0; l<ntfct; l++) { Bfct(0, k+ks, j+js, i+is) = dens0(0, k+ks, j+js, i+is)/2.;}
    
  }
};




// HOW DO WE SPECIALIZE TO SINGLE COMPONENT?
// SOME OF THESE DENSFCT BECOME DYNAMICALLY INACTIVE IE THEIR FUNCTIONAL DERIVS ARE EQUAL TO O
// THIS IS THE SAME WAY THAT CONSTANT KAPPA WORKS, ACTUALLY...

// this is rho, no-p variant
// class ce_Hs {
// public:
//   bool is_initialized;
// 
//    ce_rho_Hs() {
//      this->is_initialized = false;
// }
// 
// void YAKL_INLINE compute_P(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
// {
//   P = rho * g * z
// };
// 
// void YAKL_INLINE compute_I(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
// {
//   alpha = 1./rho;
//   theta = Theta/rho;
//   qd = (rho - \sum_n \rho_s)/rho;
//   qv = rho_v/rho;
//   ql = rho_l/rho;
//   qi = rho_i/rho;
//   I = rho * thermo.U(alpha, theta, qd, qv, ql, qi);
// };
// 
// void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
//  {
//    // ADD OR REPLACE MODE?
//    alpha = 1./rho;
//    theta = Theta/rho;
//    qd = (rho - \sum_n \rho_s)/rho;
//    qv = rho_v/rho;
//    ql = rho_l/rho;
//    qi = rho_i/rho;
//    P = rho * g * z;
//    Pi = thermo.compute_dUdtheta(alpha, theta, qd, qv, ql, qi);
//    p = -thermo.compute_dUdalpha(alpha, theta, qd, qv, ql, qi);
//    mu_d = thermo.compute_dUdqd(alpha, theta, qd, qv, ql, qi);
//    mu_v = thermo.compute_dUdqv(alpha, theta, qd, qv, ql, qi);
//    mu_l = thermo.compute_dUdql(alpha, theta, qd, qv, ql, qi);
//    mu_i = thermo.compute_dUdqi(alpha, theta, qd, qv, ql, qi);
//    dHdrho = thermo.U(alpha, theta, qd, qv, ql, qi) + P + p(alpha, theta, qd, qv, ql, qi) * alpha - theta * Pi + QN-STUFF;
//    dHdtheta = Pi;
//    dHdrhos = QN-STUFF;
//  };
// 
// };

//ADD rho_d variants
//ADD p-variants: HERE HOW TO SOLVE FOR P IS LESS CLEAR...MAYBE THIS BECOMES A FUNCTION ITSELF?

#endif
