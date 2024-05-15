#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "src/fem.h"
#include "src/mesh.h"
#include "src/solver.h"
#include "src/tests.h"
#include "src/simu.h"

/* Global variables */
std::vector< std::string > arguments;

/* To parse command line arguments */
bool flag_is_used(
    const std::string& flag,
    const std::vector< std::string >& arguments )
{
    for( int i = 0; i < arguments.size(); ++i ) {
        if( flag == arguments[i] ) {
            return true;
        }
    }
    return false;
}

using namespace FEM2A;

void run_tests()
{
    const bool t_opennl = true;
    const bool t_lmesh = true;
    const bool t_io = true;
    const bool t_Quadrature = true;
    const bool t_elem_mapping = true;
    const bool t_shape_function = true;
    const bool t_assemble_Ke = true;


    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    if( t_Quadrature ) Tests::test_Quadrature();
    if( t_elem_mapping ) Tests::test_elem_mapping();
    if( t_elem_mapping ) Tests::test_jacobian_edge();
    if( t_elem_mapping ) Tests::test_jacobian_triangle();
    if( t_shape_function ) Tests::test_shape_function();
    if( t_assemble_Ke ) Tests::test_assemble_Ke_K();
    if( t_assemble_Ke ) Tests::test_assemble_elementary_vector();

}

void run_simu()
{
    //choose the quadrature degree
    std::cout << "choose the quadrature degree (0, 2, 4, 6) : ";
    int quad_degree = 2;
    std::cin >> quad_degree;

    //choose the source function
    double (*source)(vertex) = FEM2A::Simu::sinus_bump;
    
    //test de la validité du choix de quadrature
    std::set<int> mySet {0, 2, 4, 6};
    if(mySet.find(quad_degree) == mySet.end()){
     std::cout << "!!!!!!!!!!!!!!!!     La valeur de quadrature choisie n'est pas valide ahhhhhh\n";
    }
    
    const bool simu_pure_dirichlet = false;
    const bool simu_pure_dirichlet_source = true;
    const bool simu_neumann = false;
    const bool simu_mug = false;

    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    if( simu_pure_dirichlet ) {
        Simu::pure_dirichlet_pb("data/square.mesh", verbose, quad_degree);
    }
    if( simu_pure_dirichlet_source ) {
        Simu::pure_dirichlet_pb_source("data/square.mesh", verbose, quad_degree, source);
    }
    if( simu_neumann ) {
        Simu::neumann("data/square_fine.mesh", verbose, quad_degree);
    }
    if( simu_mug ) {
        Simu::mug("data/mug_1.mesh", verbose, quad_degree);
    }
}

int main( int argc, const char * argv[] )
{
    /* Command line parsing */
    for( int i = 1; i < argc; ++i ) {
        arguments.push_back( std::string(argv[i]) );
    }

    /* Show usage if asked or no arguments */
    if( arguments.size() == 0 || flag_is_used("-h", arguments)
        || flag_is_used("--help", arguments) ) {
        std::cout << "Usage: ./fem2a [options]" << std::endl
            << "Options: " << std::endl;
        std::cout << " -h, --help:        show usage" << std::endl;
        std::cout << " -t, --run-tests:   run the tests" << std::endl;
        std::cout << " -s, --run-simu:    run the simulations" << std::endl;
        std::cout << " -v, --verbose:     print lots of details" << std::endl;
        return 0;
    }

    /* Run the tests if asked */
    if( flag_is_used("-t", arguments)
        || flag_is_used("--run-tests", arguments) ) {
        run_tests();
    }

    /* Run the simulation if asked */
    if( flag_is_used("-s", arguments)
        || flag_is_used("--run-simu", arguments) ) {
        run_simu();
    }

    return 0;
}
