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
        
        double sinus_bump( vertex v )
        {
            return 2*pow(M_PI,2)*sin(M_PI*v.x)*sin(M_PI*v.y);
        }
        
        double sin_fct( vertex v )
        {
            return sin(M_PI*v.y);
        }
        
        double sol_analytique( vertex v )
        {
            return sin(M_PI*v.x)*sin(M_PI*v.y);
        }

        //#################################
        //  Simulations
        //#################################

	// Simu 1 --------------------------------------------------------------------------------------------------

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose ){
         std::cout << "Solving a pure Dirichlet problem" << std::endl;
         int quadrature_degree = 2;
         if ( verbose ) {
          std::cout << " quadrature degree" << quadrature_degree << "\n";
          std::cout << " condition de dirichlet sur tous les bords de valeur u = x+y sans source" << "\n";
          std::cout << "maillage : " << mesh_filename << "\n";    
         }
         Mesh mesh;
         mesh.load(mesh_filename);
            
         mesh.set_attribute(unit_fct, 0, true);
            
	 Quadrature quadrature;
	 quadrature = quadrature.get_quadrature(quadrature_degree, false);
	    
	 ShapeFunctions reference_functions = ShapeFunctions(2,1);
	    
	 SparseMatrix K( mesh.nb_vertices());
	    
	 std::vector< double > F(mesh.nb_vertices(), 0);
	    
	 for (int i = 0; i < mesh.nb_triangles(); i++){
	  ElementMapping element(mesh, false, i);
	  DenseMatrix Ke;
	  Ke.set_size(3,3);
	  assemble_elementary_matrix(element, reference_functions, quadrature, unit_fct, Ke);     
       	  local_to_global_matrix(mesh, i, Ke, K);
         }
            
         std::vector< double > values(mesh.nb_vertices());
         for (int i =0 ; i < mesh.nb_vertices(); i++){
          values[i] = xy_fct(mesh.get_vertex(i)); 
         }
            
         std::vector< bool > attribute_is_dirichlet(1,true);
            
         apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
            
         std::vector<double> x(mesh.nb_vertices(), 0);
         bool converged = solve(K, F, x);
         if (!converged) {std::cout << "\nAh! le système linéaire n'a pas convergé!!!!!!!!!!!!!!!!!!!\n";}
          for (double i :x){
           std::cout << i << ' ';
          }
         std::cout <<'\n';
         //choose the name of the solution mesh saved
         save_solution(x, "square.bb");
        }

        // Simu 2 --------------------------------------------------------------------------------------------------

        void neumann_pb( const std::string& mesh_filename, bool verbose ){
         std::cout << "Solving a neumann and dirichlet problem" << std::endl;
         int quadrature_degree = 2;
	 //choose the source function (sinus_bump or unit_fct)
         auto source = unit_fct;
         auto fct = zero_fct;
         if ( verbose ) {
          std::cout << "quadrature degree : " << quadrature_degree << "\n";
          std::cout << "condition de dirichlet sur tous les bords de valeur u = " << "zero_fct" << " avec source unitaire" << "\n";
          std::cout << "maillage : " << mesh_filename << "\n";
          std::cout << "terme source : " << source << "\n";  
         }
         Mesh mesh;
         mesh.load(mesh_filename);
            
         mesh.set_attribute(unit_fct, 0, true);
            
         Quadrature quadrature;
	 quadrature = quadrature.get_quadrature(quadrature_degree, false);
	    
	 ShapeFunctions reference_functions = ShapeFunctions(2,1);
	    
	 SparseMatrix K( mesh.nb_vertices() );
	    
	 std::vector< double > F(mesh.nb_vertices() , 0);
  
	 for (int i = 0; i < mesh.nb_triangles(); i++){
	  ElementMapping element(mesh, false, i);
	      
          //K
	  DenseMatrix Ke;
          Ke.set_size(3,3);
	  assemble_elementary_matrix(element, reference_functions, quadrature, unit_fct, Ke);     
       	  local_to_global_matrix(mesh, i, Ke, K );
		
	   //F
	   std::vector< double > Fe;
	   assemble_elementary_vector(element, reference_functions, quadrature, source, Fe);     
       	   local_to_global_vector(mesh, true, i, Fe, F);
          }
            
         std::vector< double > values(mesh.nb_vertices());
         for (int i =0 ; i < mesh.nb_vertices(); i++){
          values[i] = fct(mesh.get_vertex(i)); 
         }
                
         std::vector< bool > attribute_is_dirichlet(1,true);
            
         apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
            
         std::vector<double> x(mesh.nb_vertices(), 0);
         bool converged = solve(K, F, x);
         if (!converged) {std::cout << "\nAh! le système linéaire n'a pas convergé!!!!!!!!!!!!!!!!!!!\n";}
         std::cout << "\nVecteur température :\n";
         for (double i :x){
          std::cout << i << ' ';
         }
         std::cout <<'\n';
         //choose the name of the solution mesh saved
         save_solution(x, "square.bb");
        }

        // Simu 3 --------------------------------------------------------------------------------------------------
        void pure_dirichlet_pb_source( const std::string& mesh_filename, bool verbose ){
         std::cout << "Solving a pure Dirichlet problem" << std::endl;
         int quadrature_degree = 0;
	 //choose the source function (sinus_bump or unit_fct)
         auto source = unit_fct;
         auto fct = zero_fct;
         if ( verbose ) {
          std::cout << "quadrature degree : " << quadrature_degree << "\n";
          std::cout << "condition de dirichlet sur tous les bords de valeur u = " << "zero_fct" << " avec source unitaire" << "\n";
          std::cout << "maillage : " << mesh_filename << "\n";
          std::cout << "terme source : " << source << "\n"; 
         }
         Mesh mesh;
         mesh.load(mesh_filename);
            
         mesh.set_attribute(unit_fct, 0, true);
            
         Quadrature quadrature;
	 quadrature = quadrature.get_quadrature(quadrature_degree, false);
	    
	 ShapeFunctions reference_functions = ShapeFunctions(2,1);
	    
	 SparseMatrix K( mesh.nb_vertices() );
	    
	 std::vector< double > F(mesh.nb_vertices() , 0);
	    
	 for (int i = 0; i < mesh.nb_triangles(); i++){
	  ElementMapping element(mesh, false, i);
	     
	  //K
	  DenseMatrix Ke;
	  Ke.set_size(3,3);
	  assemble_elementary_matrix(element, reference_functions, quadrature, unit_fct, Ke);     
       	  local_to_global_matrix(mesh, i, Ke, K);
		
	  //F
	  std::vector< double > Fe;
	  assemble_elementary_vector(element, reference_functions, quadrature, source, Fe);     
       	  local_to_global_vector(mesh, false, i, Fe, F);
         }
            
         std::vector< double > values(mesh.nb_vertices());
         for (int i =0 ; i < mesh.nb_vertices(); i++){
          values[i] = fct(mesh.get_vertex(i)); 
         }
             
         std::vector< bool > attribute_is_dirichlet(1,true);
            
         apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
            
         std::vector<double> x(mesh.nb_vertices(), 0);
         bool converged = solve(K, F, x);
         if (!converged) {std::cout << "\nAh! le système linéaire n'a pas convergé!!!!!!!!!!!!!!!!!!!\n";}
         std::cout << "\nVecteur température :\n";
         for (double i :x){
          std::cout << i << ' ';
         }
         std::cout <<'\n';
            
         //save p for analytical error for sinus bump and x for simulation
         std::vector<double> p(mesh.nb_vertices(), 0);
            
         for (int i; i<mesh.nb_vertices();i++){
          p[i] = x[i] - sol_analytique(mesh.get_vertex(i));
         }
            
         for (double i :p){
          std::cout << i << ' ';
         }
         std::cout <<'\n';
         //choose the name of the solution mesh saved
         save_solution(x, "square.bb");
        }
   }
}
