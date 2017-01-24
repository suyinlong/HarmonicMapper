/*!
*      \file HarmonicMapper.h
*      \brief Algorithm for harmonic mapping
*	   \author David Gu
*      \date Document 10/07/2010
*
*		Simple harmonic map, that maps a topological disk to the unit disk 
*		with minimal harmonic energy.
*/

// Simple harmonic map, that maps topological disk to a unit disk using harmonic map.

#ifndef _HARMONIC_MAPPER_H_
#define _HARMONIC_MAPPER_H_

#include "Mesh/iterators.h"
#include "Eigen/Sparse"
#include <vector>
#include "HarmonicMapperMesh.h"

#ifndef PI
#define PI 3.1415926535
#endif

namespace MeshLib
{

	template<typename M>
	class CHarmonicMapper
	{
	public:
		/*!	CHarmonicMapper constructor
		 *	\param pMesh the input mesh
		 */
		CHarmonicMapper( M* pMesh);
		/*!	CHarmonicMapper destructor
		 */
		~CHarmonicMapper();
		/*!  Compute the harmonic map using direct method
		 */
		void _map();
		/*!	Iterative method compute harmonic map
		 *	\param epsilon error threshould
		 */	
		void _iterative_map( double threshould = 1e-6 );

	protected:
		/*!	fix the boundary vertices to the unit circle
		 *  using arc length parameter
		 */
		void _set_boundary();
		/*!
		 *	compute edge weight
		 */
		void _calculate_edge_weight();

	protected:
		/*!	The input surface mesh
		 */
		M* m_pMesh;
		/*!	The boundary of m_pMesh
		 */
		typename M::CBoundary m_boundary;
		
		/*!	number of interior vertices
		 */ 
		int m_interior_vertices;
		/*! number of boundary vertices
		*/
		int m_boundary_vertices;
		/*! cosine law, return A
		 */
		double _inverse_cosine_law( double a, double b, double c );
	};

template<typename M>
double CHarmonicMapper<M>::_inverse_cosine_law( double a, double b, double c )
{
	 double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
     assert( cs <= 1.0 && cs >= -1.0 );
     return acos( cs );    
}

template<typename M>
void CHarmonicMapper<M>::_calculate_edge_weight()
{
	//compute edge length
	for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
		M::CEdge * pE = *eiter;
		M::CVertex *v1 = m_pMesh->edgeVertex1(pE);
		M::CVertex *v2 = m_pMesh->edgeVertex2(pE);
		pE->length() = (v1->point() - v2->point()).norm();
	}

	//compute corner angle
	for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
		M::CEdge * pE = *eiter;
		if (pE->boundary())
			continue;

		M::CHalfEdge * h0 = m_pMesh->edgeHalfedge(pE, 0);
		M::CHalfEdge * h1 = m_pMesh->edgeHalfedge(pE, 1);

		double h0j = static_cast<M::CEdge *>(h0->he_next()->edge())->length();
		double h0k = static_cast<M::CEdge *>(h0->he_prev()->edge())->length();

		double h1j = static_cast<M::CEdge *>(h1->he_next()->edge())->length();
		double h1k = static_cast<M::CEdge *>(h1->he_prev()->edge())->length();

		double angle0 = _inverse_cosine_law(h0j, h0k, pE->length());
		double angle1 = _inverse_cosine_law(h1j, h1k, pE->length());

		pE->weight() = 0.5 * (1/ tan(angle0) + 1/ tan(angle1));
	}
	//compute edge weight
}

/*!	CHarmonicMapper constructor 
*	Count the number of interior vertices, boundary vertices and the edge weight
*
*/
template<typename M>
CHarmonicMapper<M>::CHarmonicMapper( M* pMesh ): m_pMesh( pMesh ), m_boundary( m_pMesh )
{

	int vid  = 0; //interior vertex ID 
	int vfid = 0; //boundary vertex ID

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;

		if( pV->boundary() )
		{
			pV->idx() = vfid ++;
		}
		else
		{
			pV->idx() = vid  ++;
		}
	}

	m_interior_vertices = vid;
	m_boundary_vertices = vfid;
	
	_calculate_edge_weight();

};

//Destructor
/*!
 *	CHarmonicMapper destructor
 */
template<typename M>
CHarmonicMapper<M>::~CHarmonicMapper()
{
};


//Set the boundary vertices to the unit circle
/*!
 *	Fix the boundary using arc length parameter
 */
template<typename M>
void CHarmonicMapper<M>::_set_boundary()
{
	//get the boundary half edge loop
	std::vector<M::CLoop*> & pLs =  m_boundary.loops();
	assert( pLs.size() == 1 );
	M::CLoop * pL = pLs[0];
	std::list<M::CHalfEdge*> & pHs = pL->halfedges();
	
	//compute the total length of the boundary
	double total_length = pL->length();
	
	//parameterize the boundary using arc length parameter
	double current_length = 0;
	for (std::list<M::CHalfEdge*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter++) {
		M::CHalfEdge *pH = *hiter;

		M::CEdge *pE = m_pMesh->halfedgeEdge(pH);

		M::CVertex *v1 = m_pMesh->edgeVertex1(pE);
		M::CVertex *v2 = m_pMesh->edgeVertex2(pE);

		double length = (v1->point() - v2->point()).norm();

		double angle = 2 * PI / total_length * current_length;
		MeshLib::CPoint2 *huv = new MeshLib::CPoint2(0.5 + 0.5 * cos(angle), 0.5 + 0.5 * sin(angle));
		v1->huv() = *huv;
		current_length += length;
	}
}

//Compute the harmonic map with the boundary condition, direct method
/*!	Compute harmonic map using direct method
*/
template<typename M>
void CHarmonicMapper<M>::_map()
{
	
	//fix the boundary
	_set_boundary();
	
	std::vector<Eigen::Triplet<double> > A_coefficients;
	std::vector<Eigen::Triplet<double> > B_coefficients;

	
	//set the matrix A
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		if( pV->boundary() ) continue;
		int vid = pV->idx();

		double sw = 0;
		for( M::VertexVertexIterator witer( pV ); !witer.end(); ++ witer )
		{
			M::CVertex * pW = *witer;
			int wid = pW->idx();
			
			M::CEdge * e = m_pMesh->vertexEdge( pV, pW );
			double w = e->weight();
			if( pW->boundary() )
			{
				B_coefficients.push_back( Eigen::Triplet<double>(vid,wid,w) );
			}
			else
			{
				A_coefficients.push_back( Eigen::Triplet<double>(vid,wid, -w) );
			}
			sw += w;
		}
		A_coefficients.push_back( Eigen::Triplet<double>(vid,vid, sw ) );
	}


	Eigen::SparseMatrix<double> A( m_interior_vertices, m_interior_vertices );
	A.setZero();

	Eigen::SparseMatrix<double> B( m_interior_vertices, m_boundary_vertices );
	B.setZero();
	A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());
	B.setFromTriplets(B_coefficients.begin(), B_coefficients.end());


	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	std::cerr << "Eigen Decomposition" << std::endl;
	solver.compute(A);
	std::cerr << "Eigen Decomposition Finished" << std::endl;
	
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}


	for( int k = 0; k < 2; k ++ )
	{
		Eigen::VectorXd b(m_boundary_vertices);

		//set boundary constraints b vector
		for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
			M::CVertex * pV = *viter;
			if (!pV->boundary())
				continue;
			int id = pV->idx();
			b(id) = pV->huv()[k];
		}

		Eigen::VectorXd c(m_interior_vertices);
		c = B * b;

		Eigen::VectorXd x = solver.solve(c);
		if( solver.info() != Eigen::Success )
		{
			std::cerr << "Waring: Eigen decomposition failed" << std::endl;
		}

		//set the images of the harmonic map to interior vertices
		for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
			M::CVertex * pV = *viter;
			if (pV->boundary())
				continue;
			int id = pV->idx();
			pV->huv()[k] = x(id);
		}
	}
}

//Compute the harmonic map with the boundary condition, iterative method
/*!	Iterative method compute harmonic map
*	\param epsilon error threshould
*/
template<typename M>
void CHarmonicMapper<M>::_iterative_map( double epsilon )
{
	//fix the boundary
	_set_boundary();

	//move interior each vertex to its center of neighbors
	for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
		M::CVertex * pV = *viter;
		if (pV->boundary())
			continue;
		MeshLib::CPoint2 *huv = new MeshLib::CPoint2(0, 0);
		pV->huv() = *huv;
	}
	
	while( true )
	{
		double error = -1e+10;
		//move interior each vertex to its center of neighbors
		for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
			M::CVertex * pV = *viter;
			if (pV->boundary())
				continue;
			double total_weight = 0;
			CPoint2 huv(0, 0);
			for (M::VertexVertexIterator vviter(pV); !vviter.end(); ++vviter) {
				M::CVertex * pW = *vviter;
				M::CEdge * pE = m_pMesh->vertexEdge(pV, pW);
				double weight = pE->weight();
				total_weight += weight;
				huv = huv + pW->huv() * weight;
			}
			huv = huv / total_weight;
			double _error = (pV->huv() - huv).norm();
			error = (_error > error) ? _error : error;
			pV->huv() = huv;
		}

		if( error < epsilon ) break;
	}	
}


};
#endif

