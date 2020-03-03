
#include "fields.h"

int main(int argc, char** argv) {
  yakl::init();

// create topology
PeriodicTopology topology;
PeriodicExchange exchange;

Field field;
field.initialize();

void initialize(Topology &topo, Exchange &exch, std::string fieldName, int ndof0, int ndof1, int ndof2 = 0, int ndof3 = 0);
void waxpy(real alpha, Field &x, Field &y);
void copy(Field & f);
void clone(Field &f, std::string fieldName);
int YAKL_INLINE get_offset(int ndof);

}
