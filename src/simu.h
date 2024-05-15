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

	//condition de Dirichlet pour la première simulation
        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //terme source pour le sinus_bump
        double sinus_bump( vertex v )
        {
            return 2*pow(M_PI,2)*sin(M_PI*v.x)*sin(M_PI*v.y);
        }

        double sin_fct( vertex v )
        {
            return sin(M_PI*v.y);
        }

        //solution analytique pour le sinus bump
        double sol_analytique( vertex v )
        {
            return sin(M_PI*v.x)*sin(M_PI*v.y);
        }

        //Conditions de Neumann
	double top (vertex v)
	{
	    if (v.y >= 0.999999999){return 1.;}
	    else {return -1.;};
	}
	
	double bottom (vertex v)
	{
	    if (v.y <= 0.000000001){return 1.;}
	    else {return -1.;};
	}
	
	double right (vertex v)
	{
	    if (std::abs(v.x - 1.) <= 0.000000001){return 1.;}
	    else {return -1.;};
	}
	
	double left (vertex v)
	{
	    if (std::abs(v.x) <= 0.000000001){return 1.;}
	    else {return -1.;};
	}
	
	//condition de Neumann pour la simulation carré
	double neumann_square (vertex v)
	{
	    if (left(v) ==1){
	     return sin(M_PI*v.y);
	    } 
	    else {return 0;}
	}
	
	//conditions de Neumann et Dirichlet pour le mug
	double liquid (vertex v)
	{
	    if (((v.x == 1) && (v.y >= 1) && (v.y <= 10)) ||
		((v.x == 20) && (v.y >= 1) && (v.y <= 10)) ||
		((v.y == 1) && (v.x >= 1) && (v.x <= 20)))
		{return 1.;}
	    else {return -1.;}
		
	}
	
	double air (vertex v) //Region a lair libre
	{
	    if (((v.x == 1) && (v.y >= 1) && (v.y <= 10)) ||
		 ((v.x == 20) && (v.y >= 1) && (v.y <= 10)) ||
		 ((v.y == 1) && (v.x >= 1) && (v.x <= 20)))
            {return -1.;}
	    else {return 1.;}
		
	}
	
	double neumann_flux (vertex v)
	{
		return -0.1;
	}

        //#################################
        //  Simulations
        //#################################

	// Simu 1 --------------------------------------------------------------------------------------------------

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose, int quad_degree){
         std::cout << "Solving a pure Dirichlet problem" << std::endl;
         int quadrature_degree = quad_degree;
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

        void neumann( const std::string& mesh_filename, bool verbose, int quad_degree){
         std::cout << "Solving a neumann and dirichlet problem" << std::endl;
         //choose the quadrature degree
         int quadrature_degree = quad_degree;
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
         
         Quadrature quadrature_triangle;
	 quadrature_triangle = quadrature_triangle.get_quadrature(quadrature_degree, false);
	 Quadrature quadrature_ligne;
	 quadrature_ligne = quadrature_ligne.get_quadrature(quadrature_degree, true);
	    
	 ShapeFunctions reference_functions_triangle = ShapeFunctions(2,1);
	 ShapeFunctions reference_functions_ligne = ShapeFunctions(1,1);
	    
	 SparseMatrix K (mesh.nb_vertices());
	    
	 std::vector< double > F(mesh.nb_vertices() , 0);
  
	 for (int i = 0; i < mesh.nb_triangles(); i++){
	  ElementMapping element_triangle (mesh, false, i);
	      
          //K
	  DenseMatrix Ke;
          Ke.set_size(3,3);
	  assemble_elementary_matrix(element_triangle, reference_functions_triangle, quadrature_triangle, unit_fct, Ke);     
       	  local_to_global_matrix(mesh, i, Ke, K );
		
	   //F
	   std::vector< double > Fe(3,0.);
	   assemble_elementary_vector(element_triangle, reference_functions_triangle, quadrature_triangle, source, Fe);
       	   local_to_global_vector(mesh, false, i, Fe, F);
          }
            
         std::vector< double > values(mesh.nb_vertices());
         for (int i =0 ; i < mesh.nb_vertices(); i++){
          values[i] = fct(mesh.get_vertex(i));
         }
                
         mesh.set_attribute(right, 1, true);
	 mesh.set_attribute(left, 2, true);
	 mesh.set_attribute(top, 2, true);
	 mesh.set_attribute(bottom, 2, true);

         std::vector< bool > attribute_is_dirichlet (3, false); 
         std::vector< bool > attribute_is_neumann (3, false);
	 attribute_is_dirichlet[1] = true;
         attribute_is_neumann[2] = true;
         
         apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
         
	 for (int edge = 0; edge < mesh.nb_edges(); edge++){ 
          if (attribute_is_neumann[mesh.get_edge_attribute(edge)]){
           std::vector< double > Fe2(2,0.);
           ElementMapping element_ligne (mesh, true, edge);
	   assemble_elementary_neumann_vector(element_ligne, reference_functions_ligne, quadrature_ligne, neumann_square, Fe2);
	   local_to_global_vector(mesh, true, edge, Fe2, F);
          }
         }
            
         std::vector<double> x(mesh.nb_vertices(), 0);
         bool converged = solve(K, F, x);
         if (!converged) {std::cout << "\nAh! le système linéaire n'a pas convergé!!!!!!!!!!!!!!!!!!!\n";}
         std::cout << "\nVecteur température :\n";
         for (double i :x){
          std::cout << i << ' ';
         }
         std::cout <<'\n';
         //choose the name of the solution mesh saved
         save_solution(x, "square_fine.bb");
        }

        // Simu 3 --------------------------------------------------------------------------------------------------
        void pure_dirichlet_pb_source( const std::string& mesh_filename, bool verbose, int quad_degree, double (*source_fct)(vertex)){
         std::cout << "Solving a pure Dirichlet problem" << std::endl;
         int quadrature_degree = quad_degree;
	 //choose the source function (sinus_bump or unit_fct)
         auto source = source_fct;
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
	    
	 SparseMatrix K (mesh.nb_vertices());
	    
	 std::vector< double > F(mesh.nb_vertices() , 0);
	    
	 for (int i = 0; i < mesh.nb_triangles(); i++){
	  ElementMapping element(mesh, false, i);
	     
	  //K
	  DenseMatrix Ke;
	  Ke.set_size(3,3);
	  assemble_elementary_matrix(element, reference_functions, quadrature, unit_fct, Ke);     
       	  local_to_global_matrix(mesh, i, Ke, K);
		
	  //F
	  std::vector< double > Fe(3,0.);
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
         std::cout << "end\n";
         //choose the name of the solution mesh saved
         save_solution(x, "square.bb");
        }
        
         // Simu 4 --------------------------------------------------------------------------------------------------
        void mug( const std::string& mesh_filename, bool verbose, int quad_degree){
         std::cout << "Solving a pure Dirichlet problem" << std::endl;
         int quadrature_degree = quad_degree;
	 //choose the source function (sinus_bump or unit_fct)
         auto source = zero_fct;
         auto fct = unit_fct;
         if ( verbose ) {
          std::cout << "quadrature degree : " << quadrature_degree << "\n";
          std::cout << "condition de dirichlet sur tous les bords de valeur u = " << "zero_fct" << " avec source unitaire" << "\n";
          std::cout << "maillage : " << mesh_filename << "\n";
          std::cout << "terme source : " << source << "\n"; 
         }
         Mesh mesh;
         mesh.load(mesh_filename);
         
         Quadrature quadrature_triangle;
	 quadrature_triangle = quadrature_triangle.get_quadrature(quadrature_degree, false);
	 Quadrature quadrature_ligne;
	 quadrature_ligne = quadrature_ligne.get_quadrature(quadrature_degree, true);
	    
	 ShapeFunctions reference_functions_triangle = ShapeFunctions(2,1);
	 ShapeFunctions reference_functions_ligne = ShapeFunctions(1,1);
	    
	 SparseMatrix K (mesh.nb_vertices());
	    
	 std::vector< double > F(mesh.nb_vertices() , 0);
  
	 for (int i = 0; i < mesh.nb_triangles(); i++){
	  ElementMapping element_triangle (mesh, false, i);
	      
          //K
	  DenseMatrix Ke;
          Ke.set_size(3,3);
	  assemble_elementary_matrix(element_triangle, reference_functions_triangle, quadrature_triangle, unit_fct, Ke);     
       	  local_to_global_matrix(mesh, i, Ke, K );
		
	   //F

          }
            
         std::vector< double > values(mesh.nb_vertices());
         for (int i =0 ; i < mesh.nb_vertices(); i++){
          values[i] = fct(mesh.get_vertex(i)) * 100;
         }
                
         mesh.set_attribute(liquid, 1, true);
	 mesh.set_attribute(air, 2, true);

         std::vector< bool > attribute_is_dirichlet (3, false); 
         std::vector< bool > attribute_is_neumann (3, false);
	 attribute_is_dirichlet[1] = true;
         attribute_is_neumann[2] = true;
         
         apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
         
	 for (int edge = 0; edge < mesh.nb_edges(); edge++){ 
          if (attribute_is_neumann[mesh.get_edge_attribute(edge)]){
           std::vector< double > Fe2(2,0.);
           ElementMapping element_ligne (mesh, true, edge);
	   assemble_elementary_neumann_vector(element_ligne, reference_functions_ligne, quadrature_ligne, neumann_flux, Fe2);
	   local_to_global_vector(mesh, true, edge, Fe2, F);
          }
         }
            
         std::vector<double> x(mesh.nb_vertices(), 0);
         bool converged = solve(K, F, x);
         if (!converged) {std::cout << "\nAh! le système linéaire n'a pas convergé!!!!!!!!!!!!!!!!!!!\n";}
         std::cout << "\nVecteur température :\n";
         for (double i :x){
          std::cout << i << ' ';
         }
         std::cout <<'\n';
         //choose the name of the solution mesh saved
         save_solution(x, "mug_1.bb");
       }
   }
}
