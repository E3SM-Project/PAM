#ifndef _HAMILTONIAN_H_
#define _HAMILTONIAN_H_

#include "common.h"
#include "thermo.h"

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

// IN THE END, I THINK IT IS BEST TO ATTACH THIS INFO TO HAMILTONIANS/OTHER FUNCTIONALS THEMSELVES!
// ie knowledge of D is in HK and PVPE
// knowledge of alpha, theta, qd, rho, etc. is in CE
// then create different functionals for different progsets!
// this fits anyways with general formulation; ie a choice of "eqnset" is really a choice of Hamiltonian, which determines ndensity/ndensityfct, etc.
// THERE WILL BE SOME DUPLICATION HERE FOR SWE, BUT THAT IS FINE, MAXIMAL PERFORMANCE HERE DOESN'T MATTER!
// ie get rid of all NTRACERS_FCT stuff...

class ProgSet {
public:
  bool is_initialized;
  
  ProgSet() {
    this->is_initialized = false;
  }

  void initialize()
  {
    this->is_initialized = true;
  }
  
  virtual real YAKL_INLINE compute_D(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k) {};
  //template<uint kdof> virtual real YAKL_INLINE compute_dDdDk(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k) {};
  
  virtual real YAKL_INLINE compute_alpha(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {};
  virtual real YAKL_INLINE compute_theta(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {};
  virtual real YAKL_INLINE compute_qd(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {};
  virtual real YAKL_INLINE compute_qv(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {};
  virtual real YAKL_INLINE compute_ql(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {};
  virtual real YAKL_INLINE compute_qi(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {};
  // INITIAL CONDITION SETTING?
  // VARIOUS PARTIAL DERIVATIVES?
};  



class hvS : public ProgSet {
  real YAKL_INLINE compute_D(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
    {return dens(0, k+ks, j+js, i+is);}
};


//   alpha = 1./rho;
//   theta = Theta/rho;
//   qd = (rho - \sum_n \rho_s)/rho;
//   qv = rho_v/rho;
//   ql = rho_l/rho;
//   qi = rho_i/rho;
class rhovS : public ProgSet {   
  real YAKL_INLINE compute_D(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
    {return dens(0, k+ks, j+js, i+is);}
    
    real YAKL_INLINE compute_alpha(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
      return 1./dens0(0, k+ks, j+js, i+is);
    };
    
    real YAKL_INLINE compute_theta(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
      return dens0(1, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
    };
    real YAKL_INLINE compute_qd(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
      return (dens0(0, k+ks, j+js, i+is) - densfct0(0, k+ks, j+js, i+is) - densfct0(1, k+ks, j+js, i+is) - densfct0(2, k+ks, j+js, i+is))/dens0(0, k+ks, j+js, i+is);
    };
    real YAKL_INLINE compute_qv(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
      return densfct0(0, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
    };
    real YAKL_INLINE compute_ql(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
      return densfct0(1, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
    };
    real YAKL_INLINE compute_qi(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
      return densfct0(2, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
    };    
    
}; 






struct pvpe {
  real pv=0., pe=0.;
};

class Functional_PVPE {
public:
  ProgSet *xv;
  bool is_initialized;

   Functional_PVPE() {
     this->is_initialized = false;
}

void initialize(ProgSet &xvar)
{
  this->xv = &xvar;
  this->is_initialized = true;
}

void YAKL_INLINE compute_q0f0(realArr q0, realArr f0, const realArr v, const realArr dens, const realArr densfct, const realArr coriolis, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1> hv;
  SArray<real,4> Dv;

  // compute zeta = D2 v
  compute_D2<1>(q0, v, is, js, ks, i, j, k);

  // compute hv = R h
  Dv(0) = xv->compute_D(dens, densfct, is, js, ks, i, j, k);
  Dv(1) = xv->compute_D(dens, densfct, is, js, ks, i-1, j, k);
  Dv(2) = xv->compute_D(dens, densfct, is, js, ks, i, j-1, k);
  Dv(3) = xv->compute_D(dens, densfct, is, js, ks, i-1, j-1, k);
  R(hv, Dv);
    
  // compute q0 = zeta / hv and f0 = f / hv
    q0(0, k+ks, j+js, i+is) = q0(0, k+ks, j+js, i+is) / hv(0);
    f0(0, k+ks, j+js, i+is) = coriolis(0, k+ks, j+js, i+is) / hv(0);

}

pvpe YAKL_INLINE compute_PVPE(const realArr v, const realArr dens, const realArr densfct, const realArr coriolis, int is, int js, int ks, int i, int j, int k)
{
  pvpe vals;
  real eta, q0;
  SArray<real,1> zeta, hv;
  SArray<real,4> Dv;
  
  // compute eta = D2 v + coriolis
  compute_D2<1>(zeta, v, is, js, ks, i, j, k);
  eta = zeta(0) + coriolis(0,k+ks,j+js,i+is);

  // compute hv = R h
  Dv(0) = xv->compute_D(dens, densfct, is, js, ks, i, j, k);
  Dv(1) = xv->compute_D(dens, densfct, is, js, ks, i-1, j, k);
  Dv(2) = xv->compute_D(dens, densfct, is, js, ks, i, j-1, k);
  Dv(3) = xv->compute_D(dens, densfct, is, js, ks, i-1, j-1, k);
  R(hv, Dv);

  // compute q0 = zeta / hv
  q0 = eta / hv(0);

  vals.pv = eta;
  vals.pe = 0.5 * eta * q0;
  
  return vals;
}

};








class Hamiltonian_Hk {
  
public:
  ProgSet *xv;
  Geometry<ndims,1,1,1> *primal_geometry;
  Geometry<ndims,1,1,1> *dual_geometry;
  bool is_initialized;

   Hamiltonian_Hk() {
     this->is_initialized = false;
}

void initialize(ProgSet &xvar, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
{
  this->xv = &xvar;
  this->primal_geometry = &primal_geom;
  this->dual_geometry = &dual_geom;
  this->is_initialized = true;
}

//HERE WE ASSUME THAT DENS0 HOLDS THE TOTAL DENSITY
real YAKL_INLINE compute_KE(const realArr v, const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
{
  
  real KE, PE, IE;
  SArray<real,1> h0, h0im1, h0jm1;
  SArray<real,2> U, he;
  SArray<real,2,2> h0arr;

  //compute U = H v
  compute_H<1,diff_ord> (U, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);


  //HACK THAT RELIES ON DENSITY BEING IN 1ST SLOT OF DENSVAR...
  //CAN WE GENERALIZE?
  // Yes, could do these over NDENSITY AND NDENSITY_FCT (probably require templating on these, or that they are set in model complile costants?)
  // and then use compute_D function to get it properly...
  // Compute h0 = I h needed for phi calcs
  compute_I<1,diff_ord> (h0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
  compute_I<1,diff_ord> (h0im1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i-1, j, k);
  compute_I<1,diff_ord> (h0jm1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j-1, k);

  //compute he = phi h0
  h0arr(0,0) = h0(0);
  h0arr(1,0) = h0(0);
  h0arr(0,1) = h0im1(0);
  h0arr(1,1) = h0jm1(0);
  phi(he, h0arr);

  return 1./2. * (he(0) * ( U(0) * v(0,k+ks,j+js,i+is)) +
              + he(1) * ( U(1) * v(1,k+ks,j+js,i+is)));

}


 void YAKL_INLINE compute_dKdx(realArr F, realArr K, realArr HE, const realArr v, const realArr U, const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,2,2> D0;
  SArray<real,2> he;

  //compute he = phi * h0
  //Exploits linearity of Hodge star to compute D0 = I D = I "\sum dens" = "\sum" I dens = "\sum" dens0
  D0(0,0) = xv->compute_D(dens0, densfct0, is, js, ks, i, j, k);
  D0(0,1) = xv->compute_D(dens0, densfct0, is, js, ks, i-1, j, k);
  D0(1,0) = xv->compute_D(dens0, densfct0, is, js, ks, i, j, k);
  D0(1,1) = xv->compute_D(dens0, densfct0, is, js, ks, i, j-1, k);
  phi(he, D0);
  HE(0, k+ks, j+js, i+is) = he(0);
  HE(1, k+ks, j+js, i+is) = he(1);

  //compute F = he * U
  F(0, k+ks, j+js, i+is) = U(0, k+ks, j+js, i+is) * he(0);
  F(1, k+ks, j+js, i+is) = U(1, k+ks, j+js, i+is) * he(1);

  //compute K = 1/2 * PhiT(U,V)
  compute_phiT(K, U, v, is, js, ks, i, j, k);
  K(0, k+ks, j+js, i+is) *= 0.5;
  

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
  
  
template<uint nt, uint ntfct> class Hamiltonian_TSWE_Hs {
 public:
   ProgSet *xv;
   Geometry<ndims,1,1,1> *primal_geometry;
   Geometry<ndims,1,1,1> *dual_geometry;
   bool is_initialized;
   realArr hs0;
   
    Hamiltonian_TSWE_Hs() {
      this->is_initialized = false;
 }
 
 void initialize(ProgSet &xvar, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom) 
   //const realArr &hs)
 {
   this->xv = &xvar;
   this->primal_geometry = &primal_geom;
   this->dual_geometry = &dual_geom;
   //compute_I<1, diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   this->is_initialized = true;
 }
 
 real YAKL_INLINE compute_PE(const realArr dens, const realArr densfct, const realArr hs, int is, int js, int ks, int i, int j, int k)
 {

  //P = S * hs + 1/2 S * h + sum_nt 1/2 h * t + sum_nt 1/2 h * tfct;
  SArray<real,2> dens0;
  compute_I<2,diff_ord> (dens0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

   real PE = hs(0, k+ks, j+js, i+is) * dens0(1) + 0.5 * dens0(0) * dens(1, k+ks, j+js, i+is);
   for (int l=2; l<2+nt; l++) { PE += 0.5 * dens0(0) * dens(l, k+ks, j+js, i+is);}
   for (int l=0; l<ntfct; l++) { PE += 0.5 * dens0(0) * densfct(l, k+ks, j+js, i+is);}
   return PE;
   
 }

 real YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 {
   return 0;
 }

 //HOW SHOULD HS + OTHER CONSTANTS (GEOPOTENTIAL, ETC.) BE HANDLED IN A UNIFORM WAY
 //IDEALLY THEY ARE PART OF INITIALIZE FOR HAMILTONIAN 
 // MAYBE INITIALIZE TAKES CONSTANT VARS AS ARGUMENTS?
 // THEN WE CAN ALSO DO THINGS LIKE TAKE I HS ONCE, I GEOP ONCE, ETC.
// YES DO IT LIKE THIS...

//For now all Hs just take one extra argument: hs or geop!
//So it is ok for a bit

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






template<uint nt, uint ntfct> class Hamiltonian_SWE_Hs {
 public:
   ProgSet *xv;
   Geometry<ndims,1,1,1> *primal_geometry;
   Geometry<ndims,1,1,1> *dual_geometry;
   bool is_initialized;
   realArr hs0;
   
    Hamiltonian_SWE_Hs() {
      this->is_initialized = false;
 }
 
 void initialize(ProgSet &xvar, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom) 
   //const realArr &hs)
 {
   this->xv = &xvar;
   this->primal_geometry = &primal_geom;
   this->dual_geometry = &dual_geom;
   //compute_I<1, diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   this->is_initialized = true;
 }
 
 real YAKL_INLINE compute_PE(const realArr dens, const realArr densfct, const realArr hs, int is, int js, int ks, int i, int j, int k)
 {

  //P = g * h * hs + 1/2 g * h * h + sum_nt 1/2 h * t + sum_nt 1/2 h * tfct;
  SArray<real,1> h0;
  compute_I<1,diff_ord> (h0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

   //real PE = g * hs(0, k+ks, j+js, i+is) * h0(0) + 0.5 * g * h0(0) * dens(0, k+ks, j+js, i+is);
   //for (int l=1; l<1+nt; l++) { PE += 0.5 * dens0(0) * dens(l, k+ks, j+js, i+is);}
   //for (int l=0; l<ntfct; l++) { PE += 0.5 * dens0(0) * densfct(l, k+ks, j+js, i+is);}
   //return PE;
   
 }

 real YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 {
   return 0;
 }

  void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr hs, int is, int js, int ks, int i, int j, int k)
  {
    // dens0 is h,S,t; densfct is tfct


//ELIMINATE THIS EVENTUALLY...
    //compute I hs
    SArray<real,1> hs0;
    compute_I<1, diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    
    // Compute dHdh = g h + g hs + sum_nt 1/2 t + sum_nt 1/2 tfct
    //B(0, k+ks, j+js, i+is) = g * hs(0) + g * dens0(0, k+ks, j+js, i+is);
    //for (int l=1; l<1+nt; l++) { B(0, k+ks, j+js, i+is) += dens0(l, k+ks, j+js, i+is)/2.;}
    //for (int l=0; l<ntfct; l++) { B(0, k+ks, j+js, i+is) += densfct0(l, k+ks, j+js, i+is)/2.;}
        

    //Compute dHdt = 1/2 h
    //for (int l=1; l<1+nt; l++) { B(l, k+ks, j+js, i+is) = dens0(0, k+ks, j+js, i+is)/2.;}
    
    //Compute dHdtfct = 1/2 h
    //for (int l=0; l<ntfct; l++) { Bfct(0, k+ks, j+js, i+is) = dens0(0, k+ks, j+js, i+is)/2.;}
    
  }
};







class Hamiltonian_CE_Hs {
 public:
   ProgSet *xv;
   Geometry<ndims,1,1,1> *primal_geometry;
   Geometry<ndims,1,1,1> *dual_geometry;
   bool is_initialized;
   ThermoPotential *thermo;
   
    Hamiltonian_CE_Hs() {
      this->is_initialized = false;
 }
 
 void initialize(ProgSet &xvar, ThermoPotential &thermodynamics, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
 {
   this->xv = &xvar;
   this->thermo = &thermodynamics;
   this->primal_geometry = &primal_geom;
   this->dual_geometry = &dual_geom;
   this->is_initialized = true;
 }
 

 real YAKL_INLINE compute_PE(const realArr dens, const realArr densfct, const realArr geop, int is, int js, int ks, int i, int j, int k)
 {
   
   real rho = xv->compute_D(dens, densfct, is, js, ks, i, j, k);   
   return rho * geop(0, is, js, ks, i, j, k);
   
 }
 
 real YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 {
   real rho = xv->compute_D(dens, densfct, is, js, ks, i, j, k);

   //compute_I<1, diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);


//ISSUE- THESE SHOULD REALLY ALL BE POINT VALUES ie dens0/densfct0!
//can get them easily at i,j,k using compute_I, but then arguments are SAarray instead of realArr
//So probably need multiple versions of each function?
   real alpha = xv->compute_alpha(dens, densfct, is, js, ks, i, j, k);
   real theta = xv->compute_theta(dens, densfct, is, js, ks, i, j, k);
   real qd = xv->compute_qd(dens, densfct, is, js, ks, i, j, k);
   real ql = xv->compute_ql(dens, densfct, is, js, ks, i, j, k);
   real qi = xv->compute_qi(dens, densfct, is, js, ks, i, j, k);
   real qv = xv->compute_qv(dens, densfct, is, js, ks, i, j, k);
   return rho * thermo->compute_U(alpha, theta, qd, qv, ql, qi);
 }

  void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr geop, int is, int js, int ks, int i, int j, int k)
  {

    real alpha = xv->compute_alpha(dens0, densfct0, is, js, ks, i, j, k);
    real theta = xv->compute_theta(dens0, densfct0, is, js, ks, i, j, k);
    real qd = xv->compute_qd(dens0, densfct0, is, js, ks, i, j, k);
    real ql = xv->compute_ql(dens0, densfct0, is, js, ks, i, j, k);
    real qi = xv->compute_qi(dens0, densfct0, is, js, ks, i, j, k);
    real qv = xv->compute_qv(dens0, densfct0, is, js, ks, i, j, k);
    
    real U = thermo->compute_U(alpha, theta, qd, qv, ql, qi);
    real p = -thermo->compute_dUdalpha(alpha, theta, qd, qv, ql, qi);
    real Pi = thermo->compute_dUdtheta(alpha, theta, qd, qv, ql, qi);
    real mu_d = thermo->compute_dUdqd(alpha, theta, qd, qv, ql, qi);
    real mu_v = thermo->compute_dUdqv(alpha, theta, qd, qv, ql, qi);
    real mu_l = thermo->compute_dUdql(alpha, theta, qd, qv, ql, qi);
    real mu_i = thermo->compute_dUdqi(alpha, theta, qd, qv, ql, qi);

// HERE WE ARE ASSUMING THAT OUR PREDICTED VARS ARE RHO, RHO_S
// WHAT IS BEST WAY TO GENERALIZE?
// See above, each functional/progset is separate class...

    B(0, k+ks, j+js, i+is) = geop(is, js, ks, i, j, k) + U + p * alpha - theta * Pi + qv * (mu_d - mu_v) + ql * (mu_d - mu_l) + qi * (mu_d - mu_i);
    B(1, k+ks, j+js, i+is) = Pi;
    Bfct(0, k+ks, j+js, i+is) = mu_l - mu_d;
    Bfct(1, k+ks, j+js, i+is) = mu_v - mu_d;
    Bfct(2, k+ks, j+js, i+is) = mu_i - mu_d;
    
  }
};


//ADD p-variants: HERE HOW TO SOLVE FOR P IS LESS CLEAR...MAYBE THIS BECOMES A FUNCTION ITSELF?

#endif
