#include "Mesquite_all_headers.hpp"
using namespace Mesquite;
int main(int argc, char* argv[])
{
  MsqError err;

  char mesh_file_name[256];
  
  // command line arguments
  if (argc==1 || argc>2) {
    std::cerr << "meshfile name needed as argument.\n";
    exit(EXIT_FAILURE);
  }
  else if (argc==2) {
    std::cout << "Working with mesh file: "<< mesh_file_name << "\n";
    strcpy(mesh_file_name, argv[1]);
  } 

}
