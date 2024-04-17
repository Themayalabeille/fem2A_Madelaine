#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border )
    {
        if ( border ) {
        for (int n = 0; n < 2; n++) vertices_.push_back(M.get_edge_vertex(i, n));
        }
        else{
        for (int n = 0; n < 3; n++) vertices_.push_back(M.get_triangle_vertex(i, n));
        }
        /*/
        // afficher vecteurs de vertices_
        std::cout << "[ElementMapping] constructor for element " << i << " ";
        std::cout << "(border)"<< '\n';
        for (int v = 0; v < vertices_.size();v++){
        std::cout << "x : " << vertices_[v].x << "y : " << vertices_[v].y << "\n";
        }
        /*/
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        vertex r ;
        if ( border_ ) {
        r.x = (1-x_r.x)*vertices_[0].x + x_r.x*vertices_[1].x;
        r.y = (1-x_r.x)*vertices_[0].y + x_r.x*vertices_[1].y;
        }
        else{
        r.x = (1-x_r.x-x_r.y)*vertices_[0].x + x_r.x*vertices_[1].x + x_r.y*vertices_[2].x;
        r.y = (1-x_r.x-x_r.y)*vertices_[0].y + x_r.x*vertices_[1].y + x_r.y*vertices_[2].y;
        }
        return r ;
    }

DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        DenseMatrix J ;
        if (border_){
            J.set_size(1,2);
            J.set(0, 0, vertices_[1].x-vertices_[0].x);
            J.set(0, 1, vertices_[1].y-vertices_[0].y);
        }
        else {
            J.set_size(2,2);
            J.set(0, 0, vertices_[1].x-vertices_[0].x);
            J.set(1, 0, vertices_[1].y-vertices_[0].y);
            J.set(0, 1, vertices_[2].x-vertices_[0].x);
            J.set(1, 1, vertices_[2].y-vertices_[0].y);
        }
        return J ;
    }
    
    double ElementMapping::jacobian( vertex x_r ) const
    {
        DenseMatrix J = jacobian_matrix( x_r );
        double determinant;
        if (border_){
            determinant = pow(0.5, pow(2, J.get(0, 0)) + pow(2, J.get(0, 1)));
        }
        else {
            determinant = J.det_2x2();
        }
        return determinant ;
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        if (dim != 1 && dim != 2 || order != 1){std::cout << "alerte valeur de dim ou order impossible" << "\n";}
    }

    int ShapeFunctions::nb_functions() const
    {
        return dim_+1 ;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        if (dim_ == 2){
        switch(i){
		case 0: return 1-x_r.x-x_r.y;
		case 1: return x_r.x;
		case 2: return x_r.y;
		}
        }
        if (dim_ == 1){
        switch(i){
		case 0: return 1-x_r.x;
		case 1: return x_r.x;
		}
        }
        return 0. ; // should not be reached
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        vec2 g ;
        if (dim_ == 2){
        switch(i){
		case 0:
		g.x = -1;
		g.y = -1;
		break;
		case 1:
		g.x = 1;
		g.y = 0;
		break;
		case 2:
		g.x = 0;
		g.y = 1;
		break;
		}
        }
        if (dim_ == 1){
        switch(i){
		case 0:
		g.x = -1;
		g.y = 0;
		break;
		case 1:
		g.x = 1;
		g.y = 0;
		break;
		}
        }
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& Ke_quad,
        double (*coefficient)(vertex),
        DenseMatrix& Ke )
    {
        for (int i = 0; i < Ke.height(); i++){
         for (int j = 0; j < Ke.width(); j++){
		for (int q = 0; q < Ke_quad.nb_points(); q++){
			vertex gauss_point = Ke_quad.point(q);
			DenseMatrix jacobian = elt_mapping.jacobian_matrix(gauss_point);
			vec2 gradi = jacobian.invert_2x2().transpose().mult_2x2_2(reference_functions.evaluate_grad(i, gauss_point));
			vec2 gradj = jacobian.invert_2x2().transpose().mult_2x2_2(reference_functions.evaluate_grad(j, gauss_point));
			double scal = dot(gradi,gradj);
			Ke.add(i, j, Ke_quad.weight(q) * coefficient(elt_mapping.transform(gauss_point)) * scal * elt_mapping.jacobian(gauss_point));
		}
         }
        }
    }

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {    
        for (int i = 0; i < Ke.height(); i++){
         for (int j = 0; j < Ke.width(); j++){
          K.add(M.get_triangle_vertex_index(t,i), M.get_triangle_vertex_index(t,j), Ke.get(i,j));
        }
        }
    }

	void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        //std::cout << "compute elementary vector (source term)" << '\n';
        // TODO
        
        int imax;
        double shape_i;
        
        
        imax = reference_functions.nb_functions();
        Fe.resize(imax);
        
        for (int i =0; i<imax; ++i)
            {
            Fe[i] =0;
            for(int q = 0; q< quadrature.nb_points(); ++q)
                {
                shape_i = reference_functions.evaluate(i,quadrature.point(q));
                Fe[i]+= quadrature.weight(q)*source(elt_mapping.transform(quadrature.point(q)))* shape_i * elt_mapping.jacobian(quadrature.point(q));
                }
            }
    }


    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (neumann condition)" << '\n';
        // TODO
        double sum;
        for (int i = 0; i < reference_functions_1D.nb_functions(); i++) {
        	sum = 0;
        	for (int q = 0; q < quadrature_1D.nb_points(); q++) {
        		vertex val = quadrature_1D.point(q);
        		sum += quadrature_1D.weight(q) * reference_functions_1D.evaluate(i, val) * neumann(val) * elt_mapping_1D.jacobian(val);
        		std::cout << "i : " << i << " , q : " << q << " calcul : " << sum << "\n";
        	}        	
        	Fe.push_back(sum);
    	}
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        if (border) {
        	for (int j = 0; j < Fe.size();j++){
        	 F[M.get_edge_vertex_index(i, j)] += Fe[j];
        	}
        }
        else {
        	for (int j = 0; j < Fe.size();j++){
        	 F[M.get_triangle_vertex_index(i, j)] += Fe[j];
        	}
        }
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        int P = 10000;
        std::vector<bool> node_check;
        for (int i = 0; i < M.nb_vertices(); i++){
                 node_check.push_back(false);}
        for (int i = 0; i< M.nb_edges(); i++){
         if (attribute_is_dirichlet[M.get_edge_attribute(i)]){
           int j1 = M.get_edge_vertex_index(i,0);
           if (!node_check[j1]) {
           K.add(j1, j1, P);
           F[j1] += values[j1]*P;
           node_check[j1] = true;
           }
           int j2 = M.get_edge_vertex_index(i,1);
           if (!node_check[j2]) {
           K.add(j2, j2, P);
           F[j2] += values[j2]*P;
           node_check[j2] = true;
           }
          }
         }
    }

    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
    }

}
