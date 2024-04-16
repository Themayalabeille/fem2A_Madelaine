#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"
#include "simu.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <numeric>

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        bool test_Quadrature(){
        	std::cout << "\ntest quadrature :\n";
		for (int i = 0; i < 4; i++){
			Quadrature arr;
			arr = arr.get_quadrature(i*2, false);
			std::cout << arr.nb_points() << "\n" << arr.weight(0) << "\n";
			double sum = 0;
			for (int n = 0; n < arr.nb_points(); n++) 
			{sum += arr.weight(n);}
		    	std::cout << "somme : " << sum << "\n";
		    	if (sum != 0.5) return false;}
		return true;
        }
        
        bool test_elem_mapping(){
        	std::cout << "\ntest element_mapping :\n";
        	Mesh square;
            	square.load("data/square.mesh");
        	ElementMapping triangle4 = ElementMapping(square, false, 4);
        	vertex vec1;
        	vec1.x  = 0.2;
        	vec1.y = 0.4;
        	std::cout << "mapping x y : " << triangle4.transform(vec1).x << " " << triangle4.transform(vec1).y << "\n";
        	return true;
        }
        
        bool test_shape_function(){
                std::cout << "\ntest shape_function :\n";
        	ShapeFunctions fct1 = ShapeFunctions(2,1);
        	ShapeFunctions fct2 = ShapeFunctions(3,1);
        	std::cout << "nb fonctions :" << fct1.nb_functions() << "\n";
        	vertex vec2;
        	vec2.x  = 1;
        	vec2.y = 1;
        	int i=0;
        	std::cout << "shapefunction" << i << " : " << fct1.evaluate(i, vec2) << "\n";
        	std::cout << "shapefunction grad x y " << i << " : " << fct1.evaluate_grad(i, vec2).x << " " << fct1.evaluate_grad(0, vec2).y << "\n";
        	return true;
        	
        }
        

	bool test_jacobian_edge(){
            std::cout << "jacobian edge :\n";
            Mesh mesh;
            mesh.load("data/square.mesh");
            ElementMapping element(mesh, true, 4);
            vertex vec1;
            vec1.x  = 0.2;
            vec1.y = 0.4;
            DenseMatrix J;
            J = element.jacobian_matrix(vec1);
            std::cout << J.get(0, 0) << " , " << J.get(0, 1) << "\n";
            std::cout << "det: " << element.jacobian(vec1) << "\n";
            return true;
        }
        
        bool test_jacobian_triangle(){
            std::cout << "jacobian triangle :\n";
            Mesh mesh;
            mesh.load("data/square.mesh");
            ElementMapping element(mesh, false, 4);
            vertex vec1;
            vec1.x  = 0.2;
            vec1.y = 0.4;
            DenseMatrix J;
            J = element.jacobian_matrix(vec1);
            std::cout << J.get(0,0) << " " << J.get(1,0) << "\n" << J.get(0,1) << " " << J.get(1,1) << "\n";
            std::cout << "det : " << element.jacobian(vec1) << "\n";
            return true;
        }
        
        double unit_fct( vertex v )
        {
            return 1.;
        }
        
        bool test_assemble_Ke_K(){
		std::cout << "\ntest assemblage Ke :\n";
		Quadrature quadrature;
		quadrature = quadrature.get_quadrature(2, false);
		Mesh mesh;
		mesh.load("data/square.mesh");
		ElementMapping element(mesh, false, 4);
		ShapeFunctions reference_functions = ShapeFunctions(2,1);
		DenseMatrix Ke;
		Ke.set_size(3,3);
		assemble_elementary_matrix(
		element,
		reference_functions,
		quadrature,
		unit_fct,
		Ke);
		Ke.print();
		
		std::cout << "\ntest assemblage global K :\n";
		SparseMatrix K = SparseMatrix( mesh.nb_vertices());
		local_to_global_matrix(mesh,
		4,
		Ke,
		K );
		K.print();
		return true;
        }
        
	bool test_assemble_elementary_vector()
        {
            std::cout << "\ntest assemble Fel :" << std::endl;
            Mesh mesh;
               mesh.load("data/square.mesh");
                ElementMapping element(mesh, false, 4);
                std::vector <double > Fe;
                ShapeFunctions shape(2, 1);
                Quadrature Q;
                Q = Quadrature::get_quadrature(2);
               
                assemble_elementary_vector(element, shape, Q, Simu::unit_fct, Fe );
                for (int i = 0; i < Fe.size(); ++i){
                std::cout << Fe[i] << std::endl;
                }
               
                return true;
        }
        
        bool test_dirichlet(){
        	std::cout << "\ntest dirichlet :" << std::endl;
        	Mesh mesh;
                mesh.load("data/square.mesh");
                std::vector< double > values
                values.size() = mesh.nb_vertices()
               
        	apply_dirichlet_boundary_conditions(
        	mesh,
		const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
		const std::vector< double >& values, /* size: nb of DOFs */
		SparseMatrix& K,
		std::vector< double >& F )
        	
        	return true
        }


    }
}
