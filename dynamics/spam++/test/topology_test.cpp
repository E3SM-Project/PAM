
#include "topology.h"

int main(int argc, char** argv) {
  yakl::init();

  if (ndims == 1) {

    PeriodicTopology topology(10);
    for (int i=0;i<3;i++) {
      for (int j=0;j<3;j++) {
        if (i==0 && j==0) {continue;}
        realArr vals;
        topology.create_field("test", vals, i, j);
        // ACTUALLY ADD SOME CHECKING HERE...
    }}
  }

  if (ndims == 2) {

    PeriodicTopology topology(10,10);
    for (int i=0;i<3;i++) {
      for (int j=0;j<3;j++) {
        for (int k=0;k<3;k++) {
          if (i==0 && j==0 && k==0) {continue;}
          realArr vals;
          topology.create_field("test", vals, i, j, k);
        // ACTUALLY ADD SOME CHECKING HERE...
    }}}
  }

  if (ndims == 3) {

    PeriodicTopology topology(10,10,10);
    for (int i=0;i<3;i++) {
      for (int j=0;j<3;j++) {
        for (int k=0;k<3;k++) {
          for (int l=0;l<3;l++) {
            if (i==0 && j==0 && k==0 && l==0) {continue;}
            realArr vals;
            topology.create_field("test", vals, i, j, k, l);
            // ACTUALLY ADD SOME CHECKING HERE...
    }}}}
  }


}
