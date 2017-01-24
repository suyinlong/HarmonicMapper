/*! \file Structure.h
 *  \brief algorithms for converting one structure to another
 *  \author David Gu
 *  \date   documented on 5/6/2011
 *
 *	Algorithm for structure conversions
 */

#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#include <map>
#include <vector>
#include <complex>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

namespace MeshLib
{
/*! \brief Conversion class
*
*	Algorithm for converting from one structure to another
*/
template< typename M>
  class CStructure
  {
  public:

    /*! \brief CBaseRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
	  CStructure( M * pMesh ) { m_pMesh = pMesh; };
    /*! \brief CBaseRicciFlow destructor
	 */
	  ~CStructure(){};
	
	  /*!
	   *	Convert the metric structure to angle structure
	   */
	  void _metric_2_angle( );
	  /*!
	   *	Convert the embedding structure in R3 to metric structure
	   */
	  void _embedding_2_metric();
	  /*!
	   *	Convert angle structure to vertex curvature function
	   */
	  void _angle_2_curvature( );
	  /*!
	   *	Compute Laplace-Beltrami operator from angle structure
	   */
	  /*! Compute the edge weight
	   *  \f$ w(e) = \cot \alpha + \cot \beta \f$ where \f$\alpha,\beta\f$ are 
	   *  the two corner angles against edge \f$e\f$. \f$ w(e) = \cot \alpha \f$
	   *  If the edge is on the boundary, then 
	   */

	  void _angle_2_Laplace( );
		
	
  protected:
    /*!
     *	the input mesh
	 */
    M	  * m_pMesh;

	/*!	Euclidean Cosine law
	 *   
	 * 	\param a,b,c edge lengths
	 *  return angle C
	 */
	double _cosine_law( double a, double b, double c );

  };



//Euclidean cosine law
template<typename M>
double CStructure<M>::_cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};


//Calculate corner angle
template<typename M>
void CStructure<M>::_metric_2_angle( )
{
  
	for ( M::MeshFaceIterator fiter( m_pMesh); ! fiter.end(); fiter ++ )
  {
	  M::CFace * f = *fiter;

	  M::CHalfEdge * he[3];

      he[0] = m_pMesh->faceMostCcwHalfEdge( f ); 

	  for( int i = 0; i < 3; i ++ )
      {
        he[(i+1)%3] = m_pMesh->	faceNextCcwHalfEdge(he[i]);
      }

      double l[3];
      for(int i = 0; i < 3; i ++ )
      {
		  M::CEdge * e = m_pMesh->halfedgeEdge( he[i] );
          l[i] = e->length();
      }

      for(int i = 0; i < 3; i ++ )
      {
          he[(i+1)%3]->angle() = _cosine_law( l[(i+1)%3] ,l[(i+2)%3] ,l[i] );
      }
  }
};

template<typename M>
void CStructure<M>::_embedding_2_metric( )
{
  
	for ( M::MeshEdgeIterator eiter( m_pMesh); ! eiter.end(); eiter ++ )
	{
		M::CEdge   *  e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
		e->length() = (v1->point()-v2->point()).norm();
	}
};



//Calculate vertex curvature
template<typename M>
void CStructure<M>::_angle_2_curvature( )
{

  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
      double k  = (v->boundary() )? PI: PI * 2;
	  for( M::VertexInHalfedgeIterator vh( m_pMesh, v ); !vh.end();  ++vh )
      {
		  M::CHalfEdge * he = *vh;
          k -= he->angle();
      }
      v->k() = k;
  }

};


//Calculate vertex curvature
template<typename M>
void CStructure<M>::_angle_2_Laplace( )
{

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter ++ )
  {
	  M::CEdge * e = *eiter;
	  
	  M::CHalfEdge * he = m_pMesh->edgeHalfedge( e, 0 );
	  M::CHalfEdge * pNh = m_pMesh->faceNextCcwHalfEdge( he ); 
	 
	  double wt = cos( pNh->angle() )/sin( pNh->angle() );

	  M::CHalfEdge * sh = m_pMesh->halfedgeSym( he );
	 if( sh != NULL )
	 {
		 pNh = m_pMesh->faceNextCcwHalfEdge( sh );
		 wt +=  cos( pNh->angle() )/sin( pNh->angle() );
	 } 
	 e->weight() = wt;
	  
  }
};

};


#endif  _STRUCTURE_H_