#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            Mesh mesh;
            mesh.load(mesh_filename);
            Quadrature quadrature;
	    quadrature = quadrature.get_quadrature(2, false);
	    ShapeFunctions reference_functions = ShapeFunctions(2,1);
	    SparseMatrix K = SparseMatrix( mesh.nb_vertices());
	    for (int i = 0; i < mesh.nb_triangles(); i++){
	     ElementMapping element(mesh, false, i);
	     DenseMatrix Ke;
	     Ke.set_size(3,3);
	     assemble_elementary_matrix(
	     	 element,
	    	 reference_functions,
	    	 quadrature,
	    	 unit_fct,
	         Ke);     
       	     local_to_global_matrix(mesh,
		i,
		Ke,
		K );
            }
            for (int i = 0; i < mesh.nb_edges(); i++){
             mesh.set_attribute(unit_fct, 1, true);
            }
            std::vector< double > values;
            for (int i = 0; i < mesh.nb_vertices(); i++){
             values.push_back(xy_fct(mesh.get_edge_vertex(i,0)));
             values.push_back(xy_fct(mesh.get_edge_vertex(i,1)));
            }
            std::vector< double > F;
            for (int i = 0; i < mesh.nb_vertices(); i++){
             values.push_back(0);
            }
            std::vector< bool > attribute_is_dirichlet(1,true);
            
            apply_dirichlet_boundary_conditions(
        	mesh,
		attribute_is_dirichlet, /* size: nb of attributes */
		values, /* size: nb of DOFs */
		K,
		F );
            
            std::vector<double> x;
            bool converged = solve(
            K,
            F,
            x);
            if (!converged) {std::cout << "\nAh! le système linéaire n'a pas convergé!!!!!!!!!!!!!!!!!!!\n";}
            std::cout << x[0];

            std::cout << "TO BE IMPLEMENTED !!!" << std::endl;
        }

    }

}
