#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

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
        

    }
}
