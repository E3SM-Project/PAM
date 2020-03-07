#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_


#include "common.h"
#include "fields.h"
#include "topology.h"


template<uint ndims> class Exchange {
public:

    const Topology<ndims> *topology;

    int bufsize_x, bufsize_y, bufsize_z;
    int total_dofs;
    int ndof0, ndof1, ndof2, ndof3;

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

    realArrHost haloSendBuf_Zm_host;
    realArrHost haloRecvBuf_Zm_host;
    realArr haloSendBuf_Zm;
    realArr haloRecvBuf_Zm;
    realArrHost haloSendBuf_Zp_host;
    realArrHost haloRecvBuf_Zp_host;
    realArr haloSendBuf_Zp;
    realArr haloRecvBuf_Zp;

    bool is_initialized;


    Exchange();
    Exchange( const Exchange<ndims> &exch) = delete;
    Exchange& operator=( const Exchange<ndims> &exch) = delete;
    void printinfo();
    void initialize(const Exchange<ndims> &exch);
    void initialize(const Topology<ndims> &topo, int nd0, int nd1, int nd2, int nd3);
    void pack(const Field<ndims> &field);
    void unpack(Field<ndims> &field);
    void exchange();

};

      template<uint ndims> Exchange<ndims>::Exchange()
      {
        this->is_initialized = false;
        std::cout << "CREATED EXCHANGE\n";
      }


    template<uint ndims> void Exchange<ndims>::printinfo()
    {
      std::cout << "exchange info\n" << std::flush;
      std::cout << "bufsize_x " << this->bufsize_x << " bufsize_y " << this->bufsize_y << " bufsize_z " << this->bufsize_z << "\n" << std::flush;
      std::cout << "total_dofs " << this->total_dofs << " ndof0 " << this->ndof0 << " ndof1 " << this->ndof1 << " ndof2 " << this->ndof2 << " ndof3 " << this->ndof3 << "\n" << std::flush;
    }


template<uint ndims> void Exchange<ndims>::initialize(const Exchange<ndims> &exch)
{
 initialize(*exch.topology, exch.ndof0, exch.ndof1, exch.ndof2, exch.ndof3);
}


// ALL BROKEN FOR PARALLEL
 template<uint ndims> void Exchange<ndims>::initialize(const Topology<ndims> &topo, int nd0, int nd1, int nd2, int nd3)
 {
   this->topology = &topo;


   this->ndof0 = nd0;
   this->ndof1 = nd1;
   this->ndof2 = nd2;
   this->ndof3 = nd3;

   if (ndims == 1) { this->total_dofs = this->ndof0 + this->ndof1; }
   if (ndims == 2) { this->total_dofs = this->ndof0 + 2*this->ndof1 + this->ndof2; }
   if (ndims == 3) { this->total_dofs = this->ndof0 + 3*this->ndof1 + 3*this->ndof2 + this->ndof3; }

   if (ndims == 1)
   {
     this->bufsize_x = this->topology->halosize_x*this->total_dofs;
   }

   // This duplicates corner elements- probably okay
   if (ndims == 2)
   {
     this->bufsize_x = this->topology->halosize_x*this->total_dofs*this->topology->n_cells_y;
     this->bufsize_y = this->topology->halosize_y*this->total_dofs*this->topology->n_cells_x;
   }

   // This duplicates corner elements- probably okay
   if (ndims == 3)
   {
     this->bufsize_x = this->topology->halosize_x*this->total_dofs*this->topology->n_cells_y*this->topology->n_cells_z;
     this->bufsize_y = this->topology->halosize_y*this->total_dofs*this->topology->n_cells_x*this->topology->n_cells_z;
     this->bufsize_z = this->topology->halosize_z*this->total_dofs*this->topology->n_cells_x*this->topology->n_cells_y;
   }


   this->haloSendBuf_Xm = realArr("haloSendBuf_Xm", this->bufsize_x);
   this->haloRecvBuf_Xm = realArr("haloRecvBuf_Xm", this->bufsize_x);
   this->haloSendBuf_Xm_host = this->haloSendBuf_Xm.createHostCopy();
   this->haloRecvBuf_Xm_host = this->haloRecvBuf_Xm.createHostCopy();
   this->haloSendBuf_Xp = realArr("haloSendBuf_Xp", this->bufsize_x);
   this->haloRecvBuf_Xp = realArr("haloRecvBuf_Xp", this->bufsize_x);
   this->haloSendBuf_Xp_host = this->haloSendBuf_Xp.createHostCopy();
   this->haloRecvBuf_Xp_host = this->haloRecvBuf_Xp.createHostCopy();

   if (ndims >= 2) {
     this->haloSendBuf_Ym = realArr("haloSendBuf_Ym", this->bufsize_y);
     this->haloRecvBuf_Ym = realArr("haloRecvBuf_Ym", this->bufsize_y);
     this->haloSendBuf_Ym_host = this->haloSendBuf_Ym.createHostCopy();
     this->haloRecvBuf_Ym_host = this->haloRecvBuf_Ym.createHostCopy();
     this->haloSendBuf_Yp = realArr("haloSendBuf_Yp", this->bufsize_y);
     this->haloRecvBuf_Yp = realArr("haloRecvBuf_Yp", this->bufsize_y);
     this->haloSendBuf_Yp_host = this->haloSendBuf_Yp.createHostCopy();
     this->haloRecvBuf_Yp_host = this->haloRecvBuf_Yp.createHostCopy();
   }

   if (ndims == 3) {
     this->haloSendBuf_Zm = realArr("haloSendBuf_Zm", this->bufsize_z);
     this->haloRecvBuf_Zm = realArr("haloRecvBuf_Zm", this->bufsize_z);
     this->haloSendBuf_Zm_host = this->haloSendBuf_Zm.createHostCopy();
     this->haloRecvBuf_Zm_host = this->haloRecvBuf_Zm.createHostCopy();
     this->haloSendBuf_Zp = realArr("haloSendBuf_Zp", this->bufsize_z);
     this->haloRecvBuf_Zp = realArr("haloRecvBuf_Zp", this->bufsize_z);
     this->haloSendBuf_Zp_host = this->haloSendBuf_Zp.createHostCopy();
     this->haloRecvBuf_Zp_host = this->haloRecvBuf_Zp.createHostCopy();
   }

   this->is_initialized = true;

 }

 template<uint ndims> void Exchange<ndims>::pack(const Field<ndims> &field)
 {

   int is = this->topology->is;
   int js = this->topology->js;
   int ks = this->topology->ks;

   //pack left (x-) and  right (x+)
   yakl::parallel_for("PackLeftRight", this->bufsize_x, YAKL_LAMBDA (int iGlob) {
     int ndof, ii, k, j;
     yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_x, this->topology->n_cells_z, this->topology->n_cells_y, ndof, ii, k, j);
     this->haloSendBuf_Xm(iGlob) = field.data(ndof, k+ks, j+js, ii+is+this->topology->n_cells_x-this->topology->halosize_x);
     this->haloSendBuf_Xp(iGlob) = field.data(ndof, k+ks, j+js, ii+is);
   });

   //pack down (y-) and up (y+)
   if (ndims >=2) {
     yakl::parallel_for("PackDownUp", this->bufsize_y, YAKL_LAMBDA (int iGlob) {
       int ndof, i, k, jj;
       yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_y, this->topology->n_cells_z, this->topology->n_cells_x, ndof, jj, k, i);
       this->haloSendBuf_Ym(iGlob) = field.data(ndof, k+ks, jj+js+this->topology->n_cells_y-this->topology->halosize_y, i+is);
       this->haloSendBuf_Yp(iGlob) = field.data(ndof, k+ks, jj+js, i+is);
     });
   }

   //pack bottom (z-) and top (z+)
   if (ndims ==3) {
     yakl::parallel_for("PackBottomTop", this->bufsize_z, YAKL_LAMBDA (int iGlob) {
       int ndof, i, kk, j;
       yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_z, this->topology->n_cells_y, this->topology->n_cells_x, ndof, kk, j, i);
       this->haloSendBuf_Zm(iGlob) = field.data(ndof, kk+ks+this->topology->n_cells_z-this->topology->halosize_z, j+js, i+is);
       this->haloSendBuf_Zp(iGlob) = field.data(ndof, kk+ks, j+js, i+is);
     });
   }
 }

 template<uint ndims> void Exchange<ndims>::unpack(Field<ndims> &field)
 {

   int is = this->topology->is;
   int js = this->topology->js;
   int ks = this->topology->ks;

   //unpack left (x-) and  right (x+)
   yakl::parallel_for("UnpackLeftRight", this->bufsize_x, YAKL_LAMBDA (int iGlob) {
     int ndof, ii, k, j;
     yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_x, this->topology->n_cells_z, this->topology->n_cells_y, ndof, ii, k, j);
     field.data(ndof, k+ks, j+js, ii+is-this->topology->halosize_x) = this->haloRecvBuf_Xm(iGlob);
     field.data(ndof, k+ks, j+js, ii+is+this->topology->n_cells_x) = this->haloRecvBuf_Xp(iGlob);
   });

   //unpack down (y-) and up (y+)
   if (ndims >=2) {
     yakl::parallel_for("UnpackDownUp", this->bufsize_y, YAKL_LAMBDA (int iGlob) {
       int ndof, i, k, jj;
       yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_y, this->topology->n_cells_z, this->topology->n_cells_x, ndof, jj, k, i);
       field.data(ndof, k+ks, jj+js-this->topology->halosize_y, i+is) = this->haloRecvBuf_Ym(iGlob);
       field.data(ndof, k+ks, jj+js+this->topology->n_cells_y, i+is) = this->haloRecvBuf_Yp(iGlob);
     });
   }

   //unpack bottom (z-) and top (z+)
   if (ndims ==3) {
     yakl::parallel_for("UnpackBottomTop", this->bufsize_z, YAKL_LAMBDA (int iGlob) {
       int ndof, i, kk, j;
       yakl::unpackIndices(iGlob, this->total_dofs, this->topology->halosize_z, this->topology->n_cells_y, this->topology->n_cells_x , ndof, kk, j, i);
       field.data(ndof, kk+ks-this->topology->halosize_z, j+js, i+is) = this->haloRecvBuf_Zm(iGlob);
       field.data(ndof, kk+ks+this->topology->n_cells_z, j+js, i+is) = this->haloRecvBuf_Zp(iGlob);
     });
   }
 }



 template<uint ndims> void Exchange<ndims>::exchange()
 {
   //std:: cout << "exchanging\n";
   yakl::parallel_for( this->bufsize_x , YAKL_LAMBDA (int iGlob) {
     //std::cout << iGlob << " " << haloSendBuf_Xp(iGlob) << " " << haloSendBuf_Xm(iGlob) << "\n";
     this->haloRecvBuf_Xp(iGlob) = this->haloSendBuf_Xp(iGlob);
     this->haloRecvBuf_Xm(iGlob) = this->haloSendBuf_Xm(iGlob);
   });

   yakl::parallel_for( this->bufsize_y , YAKL_LAMBDA (int iGlob) {
     this->haloRecvBuf_Yp(iGlob) = this->haloSendBuf_Yp(iGlob);
     this->haloRecvBuf_Ym(iGlob) = this->haloSendBuf_Ym(iGlob);
   });

   yakl::parallel_for( this->bufsize_z , YAKL_LAMBDA (int iGlob) {
     this->haloRecvBuf_Zp(iGlob) = this->haloSendBuf_Zp(iGlob);
     this->haloRecvBuf_Zm(iGlob) = this->haloSendBuf_Zm(iGlob);
   });
 }

#endif
