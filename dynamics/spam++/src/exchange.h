#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_


#include "common.h"
#include "fields.h"
#include "topology.h"

// Xm SendBuf holds the cells closest to x- ie the left side -> Go into Xp RecvBuf
// Xp SendBuf holds the cells closest to x+ ie the right side -> Go into Xm RecvBuf

// Xm RecvBuf holds the halo cells closest to x- ie the left side
// Xp SendBuf holds the halo cells closest to x+ ie the right side

class Exchange {
public:

    const Topology *topology;

    int bufsize_x, bufsize_y, bufsize_xy;
    int total_dofs; //total number of dofs in a "horiz" slice ie at each layer or interface
  // basedof = 0,1,2 in 2D or 0,1 in 1D ie vertices,edges,cells or vertices,cells
  // extdof = 0,1 (0 = interfaces, 1=layers)
    int basedof, extdof, ndofs;
    int _nz, _nloop, _nloop_halo;

    realArrHost haloSendBuf_Xm_host;
    realArrHost haloRecvBuf_Xm_host;
    realArr haloSendBuf_Xm;
    realArr haloRecvBuf_Xm;
    realArrHost haloSendBuf_Xp_host;
    realArrHost haloRecvBuf_Xp_host;
    realArr haloSendBuf_Xp;
    realArr haloRecvBuf_Xp;

    realArrHost haloSendBuf_Ym_host;
    realArrHost haloRecvBuf_Ym_host;
    realArr haloSendBuf_Ym;
    realArr haloRecvBuf_Ym;
    realArrHost haloSendBuf_Yp_host;
    realArrHost haloRecvBuf_Yp_host;
    realArr haloSendBuf_Yp;
    realArr haloRecvBuf_Yp;

    realArrHost haloSendBuf_XYll_host;
    realArrHost haloRecvBuf_XYll_host;
    realArr haloSendBuf_XYll;
    realArr haloRecvBuf_XYll;
    realArrHost haloSendBuf_XYul_host;
    realArrHost haloRecvBuf_XYul_host;
    realArr haloSendBuf_XYul;
    realArr haloRecvBuf_XYul;
    realArrHost haloSendBuf_XYlr_host;
    realArrHost haloRecvBuf_XYlr_host;
    realArr haloSendBuf_XYlr;
    realArr haloRecvBuf_XYlr;
    realArrHost haloSendBuf_XYur_host;
    realArrHost haloRecvBuf_XYur_host;
    realArr haloSendBuf_XYur;
    realArr haloRecvBuf_XYur;
    
    // 4 corners in 2D
    MPI_Request sReq [4];
    MPI_Request rReq [4];

    MPI_Status  sStat[4];
    MPI_Status  rStat[4];

    bool is_initialized;


    Exchange();
    Exchange( const Exchange &exch) = delete;
    Exchange& operator=( const Exchange &exch) = delete;
    void printinfo();
    void initialize(const Exchange &exch);
    void initialize(const Topology &topo, int nd0, int nd1, int nd2);
    void pack(const Field &field);
    void unpack(Field &field);
    void exchange();
    void exchange_field(Field &field);
    void exchange_x();
    void exchange_y();
    void exchange_corners();

};

    Exchange::Exchange()
      {
        this->is_initialized = false;
        std::cout << "CREATED EXCHANGE\n";
      }


    void Exchange::printinfo()
    {
      std::cout << "exchange info\n" << std::flush;
      std::cout << "bufsize_x " << this->bufsize_x << " bufsize_y " << this->bufsize_y << "\n" << std::flush;
      std::cout << "total_dofs " << this->total_dofs << " basedof " << this->basedof << " extdof " << this->extdof << " ndofs " << this->ndofs << "\n" << std::flush;
    }

void Exchange::initialize(const Exchange &exch)
{
 initialize(*exch.topology,  exch.basedof, exch.extdof, exch.ndofs);
}


 void Exchange::initialize(const Topology &topo, int bdof, int edof, int nd)
 {
   this->topology = &topo;

   this->basedof = bdof;
   this->extdof = edof;
   this->ndofs = nd;
   
   this->total_dofs = this->ndofs ;
   if (ndims == 2 && this->basedof==1) { this->total_dofs = 2*this->ndofs; } //2 edges per cell in 2D

   if (this->extdof == 0) {
     this->_nz = this->topology->ni;
     this->_nloop = this->topology->n_cells_interfaces;
     this->_nloop_halo = this->topology->n_cells_interfaces_with_halo;
   }
   if (this->extdof == 1) {
     this->_nz = this->topology->nl;
     this->_nloop = this->topology->n_cells_layers;
     this->_nloop_halo = this->topology->n_cells_layers_with_halo;
   }
   
    this->bufsize_x = this->topology->halosize_x*this->total_dofs*this->topology->n_cells_y*this->_nz;
    this->bufsize_y = this->topology->halosize_y*this->total_dofs*this->topology->n_cells_x*this->_nz;
    this->bufsize_xy = this->topology->halosize_y*this->total_dofs*this->topology->halosize_x*this->_nz;


   this->haloSendBuf_Xm = realArr("haloSendBuf_Xm", this->bufsize_x);
   this->haloRecvBuf_Xm = realArr("haloRecvBuf_Xm", this->bufsize_x);
   this->haloSendBuf_Xm_host = this->haloSendBuf_Xm.createHostCopy();
   this->haloRecvBuf_Xm_host = this->haloRecvBuf_Xm.createHostCopy();
   this->haloSendBuf_Xp = realArr("haloSendBuf_Xp", this->bufsize_x);
   this->haloRecvBuf_Xp = realArr("haloRecvBuf_Xp", this->bufsize_x);
   this->haloSendBuf_Xp_host = this->haloSendBuf_Xp.createHostCopy();
   this->haloRecvBuf_Xp_host = this->haloRecvBuf_Xp.createHostCopy();

   if (ndims == 2) {
     this->haloSendBuf_Ym = realArr("haloSendBuf_Ym", this->bufsize_y);
     this->haloRecvBuf_Ym = realArr("haloRecvBuf_Ym", this->bufsize_y);
     this->haloSendBuf_Ym_host = this->haloSendBuf_Ym.createHostCopy();
     this->haloRecvBuf_Ym_host = this->haloRecvBuf_Ym.createHostCopy();
     this->haloSendBuf_Yp = realArr("haloSendBuf_Yp", this->bufsize_y);
     this->haloRecvBuf_Yp = realArr("haloRecvBuf_Yp", this->bufsize_y);
     this->haloSendBuf_Yp_host = this->haloSendBuf_Yp.createHostCopy();
     this->haloRecvBuf_Yp_host = this->haloRecvBuf_Yp.createHostCopy();
   }
   
     if (ndims == 2) {

     this->haloSendBuf_XYll = realArr("haloSendBuf_XYll", this->bufsize_xy);
     this->haloRecvBuf_XYll = realArr("haloRecvBuf_XYll", this->bufsize_xy);
     this->haloSendBuf_XYll_host = this->haloSendBuf_XYll.createHostCopy();
     this->haloRecvBuf_XYll_host = this->haloRecvBuf_XYll.createHostCopy();
     
     this->haloSendBuf_XYul = realArr("haloSendBuf_XYul", this->bufsize_xy);
     this->haloRecvBuf_XYul = realArr("haloRecvBuf_XYul", this->bufsize_xy);
     this->haloSendBuf_XYul_host = this->haloSendBuf_XYul.createHostCopy();
     this->haloRecvBuf_XYul_host = this->haloRecvBuf_XYul.createHostCopy();
     
     this->haloSendBuf_XYlr = realArr("haloSendBuf_XYlr", this->bufsize_xy);
     this->haloRecvBuf_XYlr = realArr("haloRecvBuf_XYlr", this->bufsize_xy);
     this->haloSendBuf_XYlr_host = this->haloSendBuf_XYlr.createHostCopy();
     this->haloRecvBuf_XYlr_host = this->haloRecvBuf_XYlr.createHostCopy();
     
     this->haloSendBuf_XYur = realArr("haloSendBuf_XYur", this->bufsize_xy);
     this->haloRecvBuf_XYur = realArr("haloRecvBuf_XYur", this->bufsize_xy);
     this->haloSendBuf_XYur_host = this->haloSendBuf_XYur.createHostCopy();
     this->haloRecvBuf_XYur_host = this->haloRecvBuf_XYur.createHostCopy();
   }

   this->is_initialized = true;

 }

 void Exchange::pack(const Field &field)
 {

   int is = this->topology->is;
   int js = this->topology->js;

// CAN WE SWAP LOOP ORDERING IN UNPACK INDICES TO USUAL K,J,I?
// YES, SHOULD BE FINE....
   //pack left (x-) and  right (x+)
   yakl::parallel_for("PackLeftRight", this->bufsize_x, YAKL_LAMBDA (int iGlob) {
     int ndof, ii, k, j;
     yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_x, this->_nz, this->topology->n_cells_y, ndof, ii, k, j);
     this->haloSendBuf_Xp(iGlob) = field.data(ndof, k, j+js, ii+is+this->topology->n_cells_x-this->topology->halosize_x);
     this->haloSendBuf_Xm(iGlob) = field.data(ndof, k, j+js, ii+is);
   });

   //pack down (y-) and up (y+)
   if (ndims ==2) {
     yakl::parallel_for("PackDownUp", this->bufsize_y, YAKL_LAMBDA (int iGlob) {
       int ndof, i, k, jj;
       yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_y, this->_nz, this->topology->n_cells_x, ndof, jj, k, i);
       this->haloSendBuf_Yp(iGlob) = field.data(ndof, k, jj+js+this->topology->n_cells_y-this->topology->halosize_y, i+is);
       this->haloSendBuf_Ym(iGlob) = field.data(ndof, k, jj+js, i+is);
     });
   }

   if (ndims ==2) {
     yakl::parallel_for("PackCorners", this->bufsize_xy, YAKL_LAMBDA (int iGlob) {
       int ndof, ii, k, jj;
       yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_y, this->topology->halosize_x, this->_nz, ndof, jj, ii, k);
       this->haloSendBuf_XYll(iGlob) = field.data(ndof, k, jj+js, ii+is);
       this->haloSendBuf_XYur(iGlob) = field.data(ndof, k, jj+js+this->topology->n_cells_y-this->topology->halosize_y, ii+is+this->topology->n_cells_x-this->topology->halosize_x);
       this->haloSendBuf_XYul(iGlob) = field.data(ndof, k, jj+js+this->topology->n_cells_y-this->topology->halosize_y, ii+is);
       this->haloSendBuf_XYlr(iGlob) = field.data(ndof, k, jj+js, ii+is+this->topology->n_cells_x-this->topology->halosize_x);
     });
   }
     

 }

 void Exchange::unpack(Field &field)
 {

   int is = this->topology->is;
   int js = this->topology->js;

   //unpack left (x-) and  right (x+)
   yakl::parallel_for("UnpackLeftRight", this->bufsize_x, YAKL_LAMBDA (int iGlob) {
     int ndof, ii, k, j;
     yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_x, this->_nz, this->topology->n_cells_y, ndof, ii, k, j);
     field.data(ndof, k, j+js, ii+is-this->topology->halosize_x) = this->haloRecvBuf_Xm(iGlob);
     field.data(ndof, k, j+js, ii+is+this->topology->n_cells_x) = this->haloRecvBuf_Xp(iGlob);
   });

   //unpack down (y-) and up (y+)
   if (ndims ==2) {
     yakl::parallel_for("UnpackDownUp", this->bufsize_y, YAKL_LAMBDA (int iGlob) {
       int ndof, i, k, jj;
       yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_y, this->_nz, this->topology->n_cells_x, ndof, jj, k, i);
       field.data(ndof, k, jj+js-this->topology->halosize_y, i+is) = this->haloRecvBuf_Ym(iGlob);
       field.data(ndof, k, jj+js+this->topology->n_cells_y, i+is) = this->haloRecvBuf_Yp(iGlob);
     });
   }

   if (ndims ==2) {
     yakl::parallel_for("UnpackCorners", this->bufsize_xy, YAKL_LAMBDA (int iGlob) {
       int ndof, ii, k, jj;
       yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_y, this->topology->halosize_x, this->_nz, ndof, jj, ii, k);
       field.data(ndof, k, jj+js-this->topology->halosize_y, ii+is-this->topology->halosize_x) = this->haloRecvBuf_XYll(iGlob);
       field.data(ndof, k, jj+js+this->topology->n_cells_y, ii+is+this->topology->n_cells_x) = this->haloRecvBuf_XYur(iGlob);
       field.data(ndof, k, jj+js-this->topology->halosize_y, ii+is+this->topology->n_cells_x) = this->haloRecvBuf_XYlr(iGlob);
       field.data(ndof, k, jj+js+this->topology->n_cells_y, ii+is-this->topology->halosize_x) = this->haloRecvBuf_XYul(iGlob);
     });
   }
   

 }



 void Exchange::exchange_x()
 {
   int ierr;

   if (this->topology->nprocx > 1) {
     yakl::fence();

     //Pre-post the receives
     ierr = MPI_Irecv( haloRecvBuf_Xm_host.data() , this->bufsize_x , REAL_MPI , this->topology->x_neigh(0) , 0 , MPI_COMM_WORLD , &this->rReq[0] );
     ierr = MPI_Irecv( haloRecvBuf_Xp_host.data() , this->bufsize_x , REAL_MPI , this->topology->x_neigh(1) , 1 , MPI_COMM_WORLD , &this->rReq[1] );

     //Copy send buffers to host
     haloSendBuf_Xm.deep_copy_to(haloSendBuf_Xm_host);
     haloSendBuf_Xp.deep_copy_to(haloSendBuf_Xp_host);
     yakl::fence();

     //Send the data
     ierr = MPI_Isend( haloSendBuf_Xm_host.data() , this->bufsize_x , REAL_MPI , this->topology->x_neigh(0) , 1 , MPI_COMM_WORLD , &this->sReq[0] );
     ierr = MPI_Isend( haloSendBuf_Xp_host.data() , this->bufsize_x , REAL_MPI , this->topology->x_neigh(1) , 0 , MPI_COMM_WORLD , &this->sReq[1] );

     //Wait for the sends and receives to finish
     ierr = MPI_Waitall(2, this->sReq, this->sStat);
     ierr = MPI_Waitall(2, this->rReq, this->rStat);

     //Copy recv buffers to device
     haloRecvBuf_Xm_host.deep_copy_to(haloRecvBuf_Xm);
     haloRecvBuf_Xp_host.deep_copy_to(haloRecvBuf_Xp);
   }

   else
   {
     yakl::parallel_for( this->bufsize_x , YAKL_LAMBDA (int iGlob) {
       this->haloRecvBuf_Xp(iGlob) = this->haloSendBuf_Xm(iGlob);
       this->haloRecvBuf_Xm(iGlob) = this->haloSendBuf_Xp(iGlob);
     });
   }


 }


 void Exchange::exchange_y()
 {
   int ierr;

   if (this->topology->nprocy > 1) {

     yakl::fence();

     //Pre-post the receives
     ierr = MPI_Irecv( haloRecvBuf_Ym_host.data() , this->bufsize_y , REAL_MPI , this->topology->y_neigh(0) , 0 , MPI_COMM_WORLD , &this->rReq[0] );
     ierr = MPI_Irecv( haloRecvBuf_Yp_host.data() , this->bufsize_y , REAL_MPI , this->topology->y_neigh(1) , 1 , MPI_COMM_WORLD , &this->rReq[1] );

     //Copy send buffers to host
     haloSendBuf_Ym.deep_copy_to(haloSendBuf_Ym_host);
     haloSendBuf_Yp.deep_copy_to(haloSendBuf_Yp_host);
     yakl::fence();

     //Send the data
     ierr = MPI_Isend( haloSendBuf_Ym_host.data() , this->bufsize_y , REAL_MPI , this->topology->y_neigh(0) , 1 , MPI_COMM_WORLD , &this->sReq[0] );
     ierr = MPI_Isend( haloSendBuf_Yp_host.data() , this->bufsize_y , REAL_MPI , this->topology->y_neigh(1) , 0 , MPI_COMM_WORLD , &this->sReq[1] );

     //Wait for the sends and receives to finish
     ierr = MPI_Waitall(2, this->sReq, this->sStat);
     ierr = MPI_Waitall(2, this->rReq, this->rStat);

     //Copy recv buffers to device
     haloRecvBuf_Ym_host.deep_copy_to(haloRecvBuf_Ym);
     haloRecvBuf_Yp_host.deep_copy_to(haloRecvBuf_Yp);
   }

   else
   {
   yakl::parallel_for( this->bufsize_y , YAKL_LAMBDA (int iGlob) {
     this->haloRecvBuf_Yp(iGlob) = this->haloSendBuf_Ym(iGlob);
     this->haloRecvBuf_Ym(iGlob) = this->haloSendBuf_Yp(iGlob);
   });
  }

 }

 void Exchange::exchange_corners()
 {
   int ierr;

   if (this->topology->nprocx > 1 || this->topology->nprocy > 1) {
     yakl::fence();
     
     //Pre-post the receives
     ierr = MPI_Irecv( haloRecvBuf_XYll_host.data() , this->bufsize_xy , REAL_MPI , this->topology->ur_neigh , 0 , MPI_COMM_WORLD , &this->rReq[0] );
     ierr = MPI_Irecv( haloRecvBuf_XYur_host.data() , this->bufsize_xy , REAL_MPI , this->topology->ll_neigh , 1 , MPI_COMM_WORLD , &this->rReq[1] );
     ierr = MPI_Irecv( haloRecvBuf_XYul_host.data() , this->bufsize_xy , REAL_MPI , this->topology->lr_neigh , 2 , MPI_COMM_WORLD , &this->rReq[2] );
     ierr = MPI_Irecv( haloRecvBuf_XYlr_host.data() , this->bufsize_xy , REAL_MPI , this->topology->ul_neigh , 3 , MPI_COMM_WORLD , &this->rReq[3] );

     //Copy send buffers to host
     haloSendBuf_XYll.deep_copy_to(haloSendBuf_XYll_host);
     haloSendBuf_XYur.deep_copy_to(haloSendBuf_XYur_host);
     haloSendBuf_XYul.deep_copy_to(haloSendBuf_XYul_host);
     haloSendBuf_XYlr.deep_copy_to(haloSendBuf_XYlr_host);
     yakl::fence();

     //Send the data
     ierr = MPI_Isend( haloSendBuf_XYll_host.data() , this->bufsize_xy , REAL_MPI , this->topology->ur_neigh , 1 , MPI_COMM_WORLD , &this->sReq[0] );
     ierr = MPI_Isend( haloSendBuf_XYur_host.data() , this->bufsize_xy , REAL_MPI , this->topology->ll_neigh , 0 , MPI_COMM_WORLD , &this->sReq[1] );
     ierr = MPI_Isend( haloSendBuf_XYul_host.data() , this->bufsize_xy , REAL_MPI , this->topology->lr_neigh , 3 , MPI_COMM_WORLD , &this->sReq[2] );
     ierr = MPI_Isend( haloSendBuf_XYlr_host.data() , this->bufsize_xy , REAL_MPI , this->topology->ul_neigh , 2 , MPI_COMM_WORLD , &this->sReq[3] );
     
     //Wait for the sends and receives to finish
     ierr = MPI_Waitall(4, this->sReq, this->sStat);
     ierr = MPI_Waitall(4, this->rReq, this->rStat);

     //Copy recv buffers to device
     haloRecvBuf_XYll_host.deep_copy_to(haloRecvBuf_XYll);
     haloRecvBuf_XYur_host.deep_copy_to(haloRecvBuf_XYur);
     haloRecvBuf_XYul_host.deep_copy_to(haloRecvBuf_XYul);
     haloRecvBuf_XYlr_host.deep_copy_to(haloRecvBuf_XYlr);
   }
   
   else
   {
 yakl::parallel_for( this->bufsize_xy , YAKL_LAMBDA (int iGlob) {
   this->haloRecvBuf_XYll(iGlob) = this->haloSendBuf_XYur(iGlob);
   this->haloRecvBuf_XYur(iGlob) = this->haloSendBuf_XYll(iGlob);
   this->haloRecvBuf_XYul(iGlob) = this->haloSendBuf_XYlr(iGlob);
   this->haloRecvBuf_XYlr(iGlob) = this->haloSendBuf_XYul(iGlob);
 });
}

}
   
// EVENTUALLY WE SHOULD BE MORE CLEVER HERE IE GROUP ALL THE SEND/RECVS IN X/Y TOGETHER, ETC.
void Exchange::exchange()
 {
   exchange_x();
   if (ndims ==2) { exchange_y(); }
   if (ndims ==2) { exchange_corners(); }
}

void Exchange::exchange_field(Field &field)
{
    this->pack(field);
    this->exchange();
    this->unpack(field);
}

#endif
