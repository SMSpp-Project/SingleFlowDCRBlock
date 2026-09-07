/*--------------------------------------------------------------------------*/
/*-------------------------- File OPTvect.h ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * A bunch of little useful template inline functions for manipulating
 * scalars and (dense and sparse) vectors, comprising:
 *
 * - Scalar functions: ABS(), sgn(), min(), max(), Swap(), CeilDiv()
 *
 * - Number-returning vector functions: Norm(), OneNorm(), INFNorm(),
 *   SumV(), MaxVecV(), MinVecV(), MaxVecI(), MinVecI(), ScalarProduct(),
 *   ScalarProduct[B[B]]()
 *
 * - Sparse/dense vector transformations: Sparsify(), SparsifyT(),
 *   SparsifyAT(), Densify(), Compact()
 *
 * - Vector operations: Vect[M]Assign[B[B]](), VectSum[B[B]](),
 *   VectSubtract[B[B]](), Vect[I]Scale[B[B]](), VectAdd(), VectDiff(),
 *   VectMult(), VectDivide(), VectXcg[B[B]]()
 *
 * - Array manipulation: Merge(), ShiftVect(), RotateVect(), ShiftRVect(),
 *   RotateRVect()
 *
 * - Searching: Match(), EqualVect(), BinSearch(), BinSearch1(),
 *   BinSearch2(), HeapIns(), HeapDel()
 *
 * \version 2.40
 *
 * \date 23 - 0 - 2003
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __OPTvect
 #define __OPTvect  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OPTtypes.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace OPTtypes_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  Scalar functions, like ABS(), min(), max(), sgn() ...               --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Scalar functions
 *  @{ */

/// absolute value of x, avoiding the (sometimes costly) abs()/fabs() calls

template<class T>
inline T ABS( const T x )
{
 return( x >= T( 0 ) ? x : -x );
 }

/*--------------------------------------------------------------------------*/
/// sign of x: 1 if x > 0, -1 if x < 0, 0 if x == 0

template<class T>
inline T sgn( const T x )
{
 return( x ? ( x > T( 0 ) ? T( 1 ) : T( -1 ) ) : T( 0 ) );
 }

/*--------------------------------------------------------------------------*/
/* max() and min() are already defined on gcc 3.x */

#ifndef __GNUC__

/// the smaller of x and y

template<class T>
inline T min( const T x , const T y )
{
 return( x <= y ? x : y );
 }

/// the larger of x and y

template<class T>
inline T max( const T x , const T y )
{
 return ( x >= y ? x : y );
 }

#else
#if __GNUC__ < 3

/// the smaller of x and y

template<class T>
inline T min( const T x , const T y )
{
 return( x <= y ? x : y );
 }

/// the larger of x and y

template<class T>
inline T max( const T x , const T y )
{
 return ( x >= y ? x : y );
 }

#endif
#endif

/*--------------------------------------------------------------------------*/
/// exchanges the values of v1 and v2

template<class T>
inline void Swap( T &v1 , T &v2 )
{
 const T temp = v1;

 v1 = v2;
 v2 = temp;
 }

/*--------------------------------------------------------------------------*/
/// the ceiling of the integer division x / y
/** Returns the ceiling of (smallest integer number not smaller than) x / y.
 * Both x and y have to be of integer types, since the `%' operation is
 * used.
 *
 * @param x   the numerator
 *
 * @param y   the denominator
 *
 * @return the smallest integer value >= x / y. */

template<class T1, class T2>
inline T1 CeilDiv( const T1 x , const T2 y )
{
 T1 temp = x / y;
 if( x % y )
  temp++;

 return( temp );
 }

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Vector functions: norms, scalar products, assignments of vectors,    --*/
/*-- multiplication by a scalar and so on. The types involved in these    --*/
/*-- operations must support the elementary arithmetical operations '+',  --*/
/*-- '-', '*' and '/', as well as the assigment '='.                      --*/
/*--                                                                      --*/
/*-- Vectors are usually just a set of n (a parameter of the function)    --*/
/*-- consecutive elements of a certain type; however, in many cases it    --*/
/*-- is useful to "restrict" them to a subset of their entries. Hence,    --*/
/*-- (read-only) vectors of indices (cIndex_Set) play a special role; for --*/
/*-- a vector g, with g{B} (where B is a vector of indices) we indicate   --*/
/*-- the "restricted" vector [ g[ B[ i ] ]. A typical reason for dealing  --*/
/*-- with "restricted" vectors is that "sparse" vectors, those having few --*/
/*-- nonzeroes w.r.t. their lenght, are very common.                      --*/
/*--                                                                      --*/
/*-- Vector of indices are meant to be "infinity-terminated", i.e., an    --*/
/*-- InINF must be found immediately after the last "valid" Index, and    --*/
/*-- sometimes they are required to be ordered, typically in increasing   --*/
/*-- sense.                                                               --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--       Number-returning functions (norms, scalar products ...)        --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Number-returning vector functions
 *  @{ */

/// returns Sum{ i = 0 .. n - 1 } g[ i ]^2 = g * g

template<class T>
inline T Norm( register const T *g , register Index n )
{
 register T t = 0;
 for( ; n-- ; )
 {
  register const T tmp = *(g++);
  t += tmp * tmp;
  }

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Norm( g{B} )

template<class T>
inline T Norm( register const T *g , register cIndex_Set B )
{
 register T t = 0;
 for( register Index h ; ( h = *(B++) ) < InINF ; )
 {
  register const T tmp = g[ h ];
  t += tmp * tmp;
  }

 return( t );
 }

/*--------------------------------------------------------------------------*/
/// returns Sum{ i = 0 .. n - 1 } ABS( g[ i ] )

template<class T>
inline T OneNorm( register const T *g , register Index n )
{
 register T t = 0;
 for( ; n-- ; )
  t += ABS( *(g++) );

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// OneNorm( g{B} )

template<class T>
inline T OneNorm( register const T *g , register cIndex_Set B )
{
 register T t = 0;
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  t += ABS( g[ h ] );

 return( t );
 }

/*--------------------------------------------------------------------------*/
/// returns Max{ i = 0 .. n - 1 } ABS( g[ i ] )

template<class T>
inline T INFNorm( register const T *g , register Index n )
{
 register T t = 0;
 for( ; n-- ; )
 {
  register const T tmp = *(g++);
  if( t < tmp )
   t = tmp;
  }

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// INFNorm( g{B} )

template<class T>
inline T INFNorm( register const T *g , register cIndex_Set B )
{
 register T t = 0;
 for( register Index h ; ( h = *(B++) ) < InINF ; )
 {
  register const T tmp = g[ h ];
  if( t < tmp )
   t = tmp;
  }

 return( t );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// returns Sum{ i = 0 .. n - 1 } g[ i ]

template<class T>
inline T SumV( register const T *g , register Index n )
{
 register T t = 0;
 for( ; n-- ; )
  t += *(g++);

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// SumV( g{B} )

template<class T>
inline T SumV( register const T *g , register cIndex_Set B )
{
 register T t = 0;
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  t += g[ h ];

 return( t );
 }

/*--------------------------------------------------------------------------*/
/// returns Max{ i = 0 .. n - 1 } g[ i ]; n *must* be > 0

template<class T>
inline T MaxVecV( register const T *g , register Index n )
{
 register T max = *g;
 for( ; --n ; )
  if( *(++g) > max )
   max = *g;

 return( max );
 }

/*--------------------------------------------------------------------------*/
/// returns Min{ i = 0 .. n - 1 } g[ i ]; n *must* be > 0

template<class T>
inline T MinVecV( register const T *g , register Index n )
{
 register T min = *g;
 for( ; --n ; )
  if( *(++g) < min )
   min = *g;

 return( min );
 }

/*--------------------------------------------------------------------------*/
/// returns the *index* of the maximum element of the n-vector v

template<class T>
inline Index MaxVecI( register const T *g , register cIndex n )
{
 register T max = *g;
 register Index maxi = 0;
 for( register Index i = maxi ; ++i < n ; )
  if( *(++g) > max )
  {
   maxi = i;
   max = *g;
   }

 return( maxi );
 }

/*--------------------------------------------------------------------------*/
/// returns the *index* of the minimum element of the n-vector v

template<class T>
inline Index MinVecI( register const T *g , register cIndex n )
{
 register T min = *g;
 register Index mini = 0;
 for( register Index i = mini ; ++i < n ; )
  if( *(++g) < max )
  {
   mini = i;
   min = *g;
   }

 return( mini );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// returns Sum{ i = 0 .. n - 1 } g1[ i ] * g2[ i ]

template<class T1, class T2>
inline T1 ScalarProduct( register const T1 *g1 , register const T2 *g2 ,
                         register Index n )
{
 register T1 t = 0;
 for( ; n-- ; )
  t += (*(g1++)) * (*(g2++));

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// ScalarProduct( g1 , g2{B} )

template<class T1, class T2>
inline T1 ScalarProduct( register const T1 *g1 , const T2 *g2 ,
                         register cIndex_Set B )
{
 register T1 t = 0;
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  t += (*(g1++)) * g2[ h ];

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// ScalarProduct( g1{B} , g2{B} ), B = intersection of B1 and B2

template<class T1, class T2>
inline T1 ScalarProduct( register const T1 *g1 , register cIndex_Set B1 ,
                         register const T2 *g2 , register cIndex_Set B2 )
{
 register T1 t = 0;
 register Index h = *B1;
 register Index k = *B2;
 if( ( h < InINF ) && ( k < InINF ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == InINF )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == InINF )
      break;
     g1++;
     }
    else
    {
     t += (*(g1++)) * (*(g2++));
     if( ( k = *(++B2) ) == InINF )
      break;
     if( ( h = *(++B1) ) == InINF )
      break;
     }

 return( t );

 }  // end( ScalarProduct( g1 , B1 , g2 , B2 )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// ScalarProduct( g1{B} , g2 )

template<class T1, class T2>
inline T1 ScalarProductB( register const T1 *g1 , const T2 *g2 ,
			  register cIndex_Set B )
{
 register T1 t = 0;
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  t += g1[ h ] * (*(g2++));

 return( t );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// ScalarProduct( g1{B} , g2{B} )

template<class T1, class T2>
inline T1 ScalarProductBB( register const T1 *g1 , const T2 *g2 ,
			   register cIndex_Set B )
{
 register T1 t = 0;
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  t += g1[ h ] * g2[ h ];

 return( t );
 }

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-- Sparse/dense vector transformations: turn "sparse" vectors into      --*/
/*-- "dense" ones and vice-versa.                                         --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Sparse/dense vector transformations
 *  @{ */

/// turns a "dense" vector into a "sparse" one, dropping the zero entries
/** Turns g from a "dense" n-vector to a "sparse" one, eliminating all items
 * that are exactly == 0. The set of nonzero items is written in B, with
 * names from Bs onwards, ordered in increasing sense; note that g itself is
 * compacted "in place" to only contain the nonzero entries.
 *
 * @param g    the "dense" n-vector to be sparsified, compacted in place to
 *             contain only its (former) nonzero entries
 *
 * @param B    the vector where the (increasing) indices of the nonzero
 *             entries of g are written; it is *not* InINF-terminated by
 *             this function
 *
 * @param n    the (original) size of g
 *
 * @param Bs   the index to be assigned to the first entry of g (default 0)
 *
 * @return a pointer to the first element of B after the last index
 *         written: this can be used for computing the number of nonzeroes
 *         in the "sparsified" vector and/or for InINF-terminating B. */

template<class T>
inline Index_Set Sparsify( register T* g , register Index_Set B ,
			   register Index n , register Index Bs = 0 )
{
 for( ; n ; n-- , g++ )
  if( *g )
   *(B++) = Bs++;
  else
   break;

 if( n )
 {
  register T* tg = g++;
  for( Bs++ ; --n ; g++ , Bs++ )
   if( *g )
   {
    *(tg++) = *g;
    *(B++) = Bs;
    }
  }

 return( B );

 }  // end( Sparsify )

/*--------------------------------------------------------------------------*/
/// as Sparsify(), but with an explicit tolerance for "being nonzero"
/** As Sparsify(), but elements are considered nonzero only if they are
 * >= eps (the idea being that all elements of g are >= 0). See Sparsify()
 * for the detailed semantics of the parameters and of the return value. */

template<class T>
inline Index_Set SparsifyT( register T* g , register Index_Set B ,
			    register Index n , register const T eps ,
			    register Index Bs = 0 )
{
 for( ; n ; n-- , g++ )
  if( *g >= eps )
   *(B++) = Bs++;
  else
   break;

 if( n )
 {
  register T* tg = g++;
  for( Bs++ ; --n ; g++ , Bs++ )
   if( *g >= eps )
   {
    *(tg++) = *g;
    *(B++) = Bs;
    }
  }

 return( B );

 }  // end( SparsifyT )

/*--------------------------------------------------------------------------*/
/// as SparsifyT(), but the tolerance is compared against the absolute value
/** As SparsifyT(), but elements are considered nonzero only if their
 * ABS() is >= eps. See Sparsify() for the detailed semantics of the
 * parameters and of the return value. */

template<class T>
inline Index_Set SparsifyAT( register T* g , register Index_Set B ,
			     register Index n , register const T eps ,
			     register Index Bs = 0 )
{
 for( ; n ; n-- , g++ )
  if( ABS( *g ) >= eps )
   *(B++) = Bs++;
  else
   break;

 if( n )
 {
  register T* tg = g++;
  for( Bs++ ; --n ; g++ , Bs++ )
   if( ABS( *g ) >= eps )
   {
    *(tg++) = *g;
    *(B++) = Bs;
    }
  }

 return( B );

 }  // end( SparsifyAT )

/*--------------------------------------------------------------------------*/
/// turns a "sparse" vector into a "dense" one, padding with zeroes
/** Turns g from a "sparse" m-vector, whose set of nonzero elements is B, to
 * a "dense" n-vector padded with zeroes where necessary. B has to be
 * ordered in increasing sense, but does not need to be InINF-terminated
 * (m gives the same information). Note that the function will write in
 * g[ n - 1 ], hence the vector has to have been properly allocated.
 *
 * @param g    the vector, containing the m nonzero entries "compacted" at
 *             its beginning, to be expanded in place to a "dense" n-vector
 *
 * @param B    the (increasing) indices, in 0 .. n - 1, of the m nonzero
 *             entries of g
 *
 * @param m    the number of nonzero entries in g (== the size of B)
 *
 * @param n    the size of the "dense" vector to be produced
 *
 * @param k    if k > 0, the function only works on the subvector of g
 *             between k and n - 1, i.e., g[ 0 ] .. g[ k - 1 ] are left
 *             intact while the rest is densified; it is required that B[]
 *             only contains indices >= k (and, of course, < n). */

template<class T>
inline void Densify( register T* g , register cIndex_Set B ,
		     register Index m , register Index n , cIndex k = 0 )
{
 if( m )  // there is at least a nonzero element
 {
  B += m;
  for( register Index h = *(--B) ; n > m ; )
   if( h == --n )
   {
    g[ n ] = g[ --m ];
    if( m )
     h = *(--B);
    else
     break;
    }
   else
    g[ n ] = 0;
  }

 if( ( ! m ) && ( n > k ) )
  VectAssign( g , T( 0 ) , n - k );

 } // end( Densify )

/*--------------------------------------------------------------------------*/
/// removes from a "dense" vector all the entries whose index is in B
/** Takes a "dense" n-vector g and "compacts" it, deleting the elements
 * whose indices are in B; the remaining entries in g[] are shifted left of
 * the minimum possible amount in order to fill the holes left by the
 * deleted ones.
 *
 * @param g   the "dense" n-vector to be compacted in place
 *
 * @param B   the (increasing), InINF-terminated indices, all in the range
 *            0 .. n - 1, of the entries of g to be deleted
 *
 * @param n   the (original) size of g. */

template<class T>
inline void Compact( register T* g , register cIndex_Set B ,
		     register Index n )
{
 register Index i = *(B++);  // current position where to write
 register Index j = i + 1;   // element to copy

 for( register Index h = *(B++) ; h < InINF ; j++ )
 {
  while( j < h )
   g[ i++ ] = g[ j++ ];

  h = *(B++);
  }

 VectAssign( g + i , g + j , n - j );

 }  // end( Compact )

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-- Vector operations: sum/difference of vectors, multiplication by a    --*/
/*-- scalar, assignment, setting, scaling ...                             --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Vector operations
 *  @{ */

/// g[ i ] = x for each i = 0 .. n - 1

template<class T>
inline void VectAssign( T *const g , register const T x , cIndex n )
{
 for( register T *tg = g + n ; tg > g ; )
  *(--tg) = x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g{B} = x, all other entries of g unchanged

template<class T>
inline void VectAssign( register T *const g , register const T x ,
			register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g[ h ] = x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 := g2

template<class T1, class T2>
inline void VectAssign( register T1 *g1 , register const T2 *g2 ,
                        register Index n )
{
 for( ; n-- ; )
  *(g1++) = *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 := g2, special version for g2 an InINF-terminated vector of indices
/** Special version of VectAssign( g1 , g2 , n ) for the case where g2 is an
 * InINF-terminated vector of indices: the terminating InINF is *not*
 * written into g1.
 *
 * @return a pointer to the position in g1 right after the last index
 *         copied, where the terminating InINF should be written by the
 *         caller if needed. */

inline Index_Set VectAssign( register Index_Set g1 , register cIndex_Set g2 )
{
 for( register Index h ; ( h = *(g2++) ) < InINF ; )
  *(g1++) = h;

 return( g1 );
 }  

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 = g2{B}

template<class T1, class T2>
inline void VectAssign( register T1 *g1 , register const T2 *g2 ,
                        register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g1++) = g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = g2{B}, B = intersection of B1 and B2, all other unchanged

template<class T1, class T2>
inline void VectAssign( register T1 *g1 , register cIndex_Set B1 ,
			register const T2 *g2 , register cIndex_Set B2 )
{
 register Index h = *B1;
 register Index k = *B2;
 if( ( h < InINF ) && ( k < InINF ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == InINF )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == InINF )
      break;
     g1++;
     }
    else
    {
     (*(g1++)) = (*(g2++));
     if( ( k = *(++B2) ) == InINF )
      break;
     if( ( h = *(++B1) ) == InINF )
      break;
     }

 }  // end( VectAssign( g1 , B1 , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = g2{B}, B = intersection of B1 and B2, g1{B1 / B} = gg

template<class T1, class T2>
inline void VectAssign( register T1 *g1 , register cIndex_Set B1 ,
			register const T2 *g2 , register cIndex_Set B2 ,
			register const T1 gg )
{
 register Index k = *B2;
 for( register Index h ; ( h = *(B1++) ) < InINF ; )
 {
  while( k < h )
  {
   k = *(++B2);
   g2++;
   }

  if( k == InINF )
   break;

  if( h == k )
   *(g1++) = *(g2++);
  else
   *(g1++) = gg;
  }
 }  // end( VectAssign( g1 , B1 , g2 , B2 , gg )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 := x * g2

template<class T1, class T2, class T3>
inline void VectAssign( register T1 *g1 , register  const T2 *g2 ,
                        register const T3 x , register Index n )
{
 for( ; n-- ; )
  *(g1++) = x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 = x * g2{B}

template<class T1, class T2, class T3>
inline void VectAssign( register T1 *g1 , register const T2 *g2 ,
                        register const T3 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g1++) = x * g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = x * g2{B}, B = intersection of B1 and B2, all other unchanged

template<class T1, class T2, class T3>
inline void VectAssign( register T1 *g1 , register cIndex_Set B1 ,
			register const T3 x , register const T2 *g2 ,
			register cIndex_Set B2 )
{
 register Index h = *B1;
 register Index k = *B2;
 if( ( h < InINF ) && ( k < InINF ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == InINF )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == InINF )
      break;
     g1++;
     }
    else
    {
     (*(g1++)) = x * (*(g2++));
     if( ( k = *(++B2) ) == InINF )
      break;
     if( ( h = *(++B1) ) == InINF )
      break;
     }

 }  // end( VectAssign( g1 , B1 , x , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = x * g2{B}, B = intersection of B1 and B2, g1{B1 / B} = gg

template<class T1, class T2, class T3>
inline void VectAssign( register T1 *g1 , register cIndex_Set B1 ,
			register const T3 x , register const T2 *g2 ,
			register cIndex_Set B2 , register const T1 gg )
{
 register Index k = *B2;
 for( register Index h ; ( h = *(B1++) ) < InINF ; )
 {
  while( k < h )
  {
   k = *(++B2);
   g2++;
   }

  if( k == InINF )
   break;

  if( h == k )
   *(g1++) = x * (*(g2++));
  else
   *(g1++) = gg;
  }
 }  // end( VectAssign( g1 , B1 , x , g2 , B2 , gg )

/*--------------------------------------------------------------------------*/
/// g{B} = x, all other entries of g unchanged

template<class T>
inline void VectAssignB( register T *const g , register const T x ,
			 register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g[ h ] = x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g{B} = x, all other entries (0 .. n - 1) of g1 = gg

template<class T>
inline void VectAssignB( register T *g , register const T x ,
			 register cIndex_Set B , register cIndex n ,
			 register const T gg = 0 )
{
 register Index h = *(B++);
 for( register Index i = 0 ; i < n ; )
  if( h == i++ )
  {
   *(g++) = x;
   h = *(B++);
   }
  else
   *(g++) = gg;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = g2, all other entries of g1 unchanged

template<class T1, class T2>
inline void VectAssignB( register T1 *g1 , register const T2 *g2 ,
			 register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] = *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = g2, all other entries (0 .. n - 1) of g1 = gg

template<class T1, class T2>
inline void VectAssignB( register T1 *g1 , register const T2 *g2 ,
			 register cIndex_Set B , register cIndex n ,
			 register const T1 gg = 0 )
{
 register Index h = *(B++);
 for( register Index i = 0 ; i < n ; )
  if( h == i++ )
  {
   *(g1++) = *(g2++);
   h = *(B++);
   }
  else
   *(g1++) = gg;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = x * g2, all other entries of g1 unchanged

template<class T1, class T2>
inline void VectAssignB( register T1 *g1 , register const T2 *g2 ,
			 register const T1 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] = x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = x * g2, all other entries (0 .. n - 1) of g1 = gg

template<class T1, class T2>
inline void VectAssignB( register T1 *g1 , register const T2 *g2 ,
			 register const T1 x , register cIndex_Set B ,
			 register cIndex n , register const T1 gg = 0 )
{
 register Index h = *(B++);
 for( register Index i = 0 ; i < n ; )
  if( h == i++ )
  {
   *(g1++) = x * (*(g2++));
   h = *(B++);
   }
  else
   *(g1++) = gg;
 }

/*--------------------------------------------------------------------------*/
/// g1{B} = g2{B}, all other entries unchanged

template<class T1, class T2>
inline void VectAssignBB( register T1 *g1 , register const T2 *g2 ,
			  register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] = g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = x * g2{B}, all other entries unchanged

template<class T1, class T2>
inline void VectAssignBB( register T1 *g1 , register const T2 *g2 ,
			  register const T1 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] = x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// g1 := - g2

template<class T1, class T2>
inline void VectMAssign( register T1 *g1 , register const T2 *g2 ,
			 register Index n )
{
 for( ; n-- ; )
  *(g1++) = - *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 = - g2{B}

template<class T1, class T2>
inline void VectMAssign( register T1 *g1 , register const T2 *g2 ,
			 register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g1++) = - g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = - g2{B}, B = intersection of B1 and B2, all other unchanged

template<class T1, class T2>
inline void VectMAssign( register T1 *g1 , register cIndex_Set B1 ,
			 register const T2 *g2 , register cIndex_Set B2 )
{
 register Index h = *B1;
 register Index k = *B2;
 if( ( h < InINF ) && ( k < InINF ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == InINF )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == InINF )
      break;
     g1++;
     }
    else
    {
     (*(g1++)) = - (*(g2++));
     if( ( k = *(++B2) ) == InINF )
      break;
     if( ( h = *(++B1) ) == InINF )
      break;
     }

 }  // end( VectMAssign( g1 , B1 , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = - g2{B}, B = intersection of B1 and B2, g1{B1 / B} = gg

template<class T1, class T2>
inline void VectMAssign( register T1 *g1 , register cIndex_Set B1 ,
			 register const T2 *g2 , register cIndex_Set B2 ,
			 register const T1 gg )
{
 register Index k = *B2;
 for( register Index h ; ( h = *(B1++) ) < InINF ; )
 {
  while( k < h )
  {
   k = *(++B2);
   g2++;
   }

  if( k == InINF )
   break;

  if( h == k )
   *(g1++) = - *(g2++);
  else
   *(g1++) = gg;
  }
 }  // end( VectMAssign( g1 , B1 , g2 , B2 , gg )

/*--------------------------------------------------------------------------*/
/// g1{B} = - g2, all other entries of g1 = unchanged

template<class T1, class T2>
inline void VectMAssignB( register T1 *g1 , register const T2 *g2 ,
			  register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g1++) = - g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} = - g2, all other entries (0 .. n - 1) of g1 = gg

template<class T1, class T2>
inline void VectMAssignB( register T1 *g1 , register const T2 *g2 ,
			  register cIndex_Set B , register cIndex n ,
			  register const T1 gg = 0 )
{
 register Index h = *(B++);
 for( register Index i = 0 ; i < n ; )
  if( h == i++ )
  {
   *(g1++) = - *(g2++);
   h = *(B++);
   }
  else
   *(g1++) = gg;
 }

/*--------------------------------------------------------------------------*/
/// g1{B} = - g2{B}, all other entries unchanged

template<class T1, class T2>
inline void VectMAssignBB( register T1 *g1 , register const T2 *g2 ,
			   register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] = - g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// g[ i ] += x for each i = 0 .. n - 1

template<class T>
inline void VectSum( T *const g , register const T x , cIndex n )
{
 for( register T *tg = g + n ; tg > g ; )
  *(--tg) += x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// special version of the above for g an InINF-terminated vector of indices

inline void VectSum( register Index_Set g , cIndex x )
{
 while( *g < InINF )
  *(g++) += x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g{B} += x, all other entries of g unchanged

template<class T>
inline void VectSum( register T *const g , register const T x ,
		     register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g[ h ] += x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 += g2

template<class T1, class T2>
inline void VectSum( register T1 *g1 , register const T2 *g2 ,
                     register Index n )
{
 for( ; n-- ; )
  *(g1++) += *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 += g2{B}

template<class T1, class T2>
inline void VectSum( register T1 *g1 , register const T2 *g2 , 
                     register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g1++) += g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} += g2{B}, B = intersection of B1 and B2, all other unchanged

template<class T1, class T2>
inline void VectSum( register T1 *g1 , register cIndex_Set B1 ,
		     register const T2 *g2 , register cIndex_Set B2 )
{
 register Index h = *B1;
 register Index k = *B2;
 if( ( h < InINF ) && ( k < InINF ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == InINF )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == InINF )
      break;
     g1++;
     }
    else
    {
     (*(g1++)) += (*(g2++));
     if( ( k = *(++B2) ) == InINF )
      break;
     if( ( h = *(++B1) ) == InINF )
      break;
     }

 }  // end( VectSum( g1 , B1 , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 += x * g2

template<class T1, class T2, class T3>
inline void VectSum( register T1 *g1 , register const T2 *g2 , 
                     register const T3 x , register Index n )
{
 for( ; n-- ; )
  *(g1++) += x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 += x * g2{B}

template<class T1, class T2, class T3>
inline void VectSum( register T1 *g1 , register const T2 *g2 , 
                     register const T3 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g1++) += g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} += x * g2{B}, B = intersection of B1 and B2, all other unchanged

template<class T1, class T2, class T3>
inline void VectSum( register T1 *g1 , register cIndex_Set B1 ,
		     register const T3 x , register const T2 *g2 ,
		     register cIndex_Set B2 )
{
 register Index h = *B1;
 register Index k = *B2;
 if( ( h < InINF ) && ( k < InINF ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == InINF )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == InINF )
      break;
     g1++;
     }
    else
    {
     (*(g1++)) += x * (*(g2++));
     if( ( k = *(++B2) ) == InINF )
      break;
     if( ( h = *(++B1) ) == InINF )
      break;
     }

 }  // end( VectSum( g1 , B1 , x , g2 , B2 )

/*--------------------------------------------------------------------------*/
/// g{B} += x, all other entries (0 .. n - 1) of g1 += gg

template<class T>
inline void VectSumB( register T *g , register const T x ,
		      register cIndex_Set B , register cIndex n ,
		      register const T gg )
{
 register Index h = *(B++);
 for( register Index i = 0 ; i < n ; )
  if( h == i++ )
  {
   *(g++) += x;
   h = *(B++);
   }
  else
   *(g++) += gg;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} += g2, all other entries of g1 unchanged

template<class T1, class T2>
inline void VectSumB( register T1 *g1 , register const T2 *g2 ,
		      register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] += *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} += g2, all other entries (0 .. n - 1) of g1 += gg

template<class T1, class T2>
inline void VectSumB( register T1 *g1 , register const T2 *g2 ,
		      register cIndex_Set B , register cIndex n ,
		      register const T1 gg )
{
 register Index h = *(B++);
 for( register Index i = 0 ; i < n ; )
  if( h == i++ )
  {
   *(g1++) += *(g2++);
   h = *(B++);
   }
  else
   *(g1++) += gg;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} += x * g2, all other entries of g1 unchanged

template<class T1, class T2, class T3>
inline void VectSumB( register T1 *g1 , register const T2 *g2 ,
		      register const T3 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] += x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} += x * g2, all other entries (0 .. n - 1) of g1 += gg

template<class T1, class T2, class T3>
inline void VectSumB( register T1 *g1 , register const T2 *g2 ,
		      register const T3 x , register cIndex_Set B ,
		      register cIndex n , register const T1 gg )
{
 register Index h = *(B++);
 for( register Index i = 0 ; i < n ; )
  if( h == i++ )
  {
   *(g1++) += x * (*(g2++));
   h = *(B++);
   }
  else
   *(g1++) += gg;
 }

/*--------------------------------------------------------------------------*/
/// g1{B} += g2{B}, all other entries unchanged

template<class T1, class T2>
inline void VectSumBB( register T1 *g1 , register const T2 *g2 ,
		       register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] += g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} += x * g2{B}, all other entries unchanged

template<class T1, class T2, class T3>
inline void VectSumBB( register T1 *g1 , register const T2 *g2 ,
		       register const T3 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] += x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// g[ i ] -= x for each i = 0 .. n - 1
/** Subtracts the scalar x from each entry of g. Useful for *unsigned* data
 * types, for which VectSum( g , -x , n ) would not work. */

template<class T>
inline void VectSubtract( T *const g , register const T x , cIndex n )
{
 for( register T *tg = g + n ; tg > g ; )
  *(--tg) -= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// as VectSubtract( g , x , n ), for g an InINF-terminated vector of indices
/** Special version of VectSubtract( g , x , n ) for the case where g is an
 * InINF-terminated vector of indices (note that indices are generally
 * unsigned, hence the same caveat about VectSum( g , -x , n ) applies). */

inline void VectSubtract( register Index_Set g , cIndex x )
{
 while( *g < InINF )
  *(g++) -= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g{B} -= x, all other entries of g unchanged

template<class T>
inline void VectSubtract( T *const g , register const T x ,
			  register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g[ h ] -= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 -= g2 (element-wise)

template<class T1, class T2>
inline void VectSubtract( register T1 *g1 , register const T2 *g2 ,
                          register Index n )
{
 for( ; n-- ; )
  *(g1++) -= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 -= g2{B} (element-wise)

template<class T1, class T2>
inline void VectSubtract( register T1 *g1 , register const T2 *g2 , 
                          register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g1++) -= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} -= g2{B}, B = intersection of B1 and B2, all other unchanged

template<class T1, class T2>
inline void VectSubtract( register T1 *g1 , register cIndex_Set B1 ,
			  register const T2 *g2 , register cIndex_Set B2 )
{
 register Index h = *B1;
 register Index k = *B2;
 if( ( h < InINF ) && ( k < InINF ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == InINF )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == InINF )
      break;
     g1++;
     }
    else
    {
     (*(g1++)) -= (*(g2++));
     if( ( k = *(++B2) ) == InINF )
      break;
     if( ( h = *(++B1) ) == InINF )
      break;
     }

 }  // end( VectSubtract( g1 , B1 , g2 , B2 )

/*--------------------------------------------------------------------------*/
/// g{B} -= x, all other entries (0 .. n - 1) of g1 -= gg

template<class T>
inline void VectSubtractB( register T *g , register const T x ,
			   register cIndex_Set B , register cIndex n ,
			   register const T gg )
{
 register Index h = *(B++);
 for( register Index i = 0 ; i < n ; )
  if( h == i++ )
  {
   *(g++) -= x;
   h = *(B++);
   }
  else
   *(g++) -= gg;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} -= g2, all other entries of g1 unchanged

template<class T1, class T2>
inline void VectSubtractB( register T1 *g1 , register const T2 *g2 ,
			   register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] -= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} -= g2, all other entries (0 .. n - 1) of g1 -= gg

template<class T1, class T2>
inline void VectSubtractB( register T1 *g1 , register const T2 *g2 ,
			   register cIndex_Set B , register cIndex n ,
			   register const T1 gg )
{
 register Index h = *(B++);
 for( register Index i = 0 ; i < n ; )
  if( h == i++ )
  {
   *(g1++) -= *(g2++);
   h = *(B++);
   }
  else
   *(g1++) -= gg;
 }

/*--------------------------------------------------------------------------*/
/// g1{B} -= g2{B}, all other entries unchanged

template<class T1, class T2>
inline void VectSubtractBB( register T1 *g1 , register const T2 *g2 ,
			    register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] -= g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// g *= x, x a scalar

template<class T>
inline void VectScale( register T *g , register const T x , register Index n )
{
 for( ; n-- ; )
  *(g++) *= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g{B} *= x, x a scalar

template<class T>
inline void VectScale( register T *g , register const T x ,
		       register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g[ h ] *= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 *= g2

template<class T1, class T2>
inline void VectScale( register T1 *g1 , register const T2 *g2 ,
		       register Index n )
{
 for( ; n-- ; )
  *(g1++) *= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 *= g2{B}

template<class T1, class T2>
inline void VectScale( register T1 *g1 , register const T2 *g2 , 
		       register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g1++) *= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} *= g2{B}, B = intersection of B1 and B2, all other unchanged

template<class T1, class T2>
inline void VectScale( register T1 *g1 , register cIndex_Set B1 ,
		       register const T2 *g2 , register cIndex_Set B2 )
{
 register Index h = *B1;
 register Index k = *B2;
 if( ( h < InINF ) && ( k < InINF ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == InINF )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == InINF )
      break;
     g1++;
     }
    else
    {
     (*(g1++)) *= (*(g2++));
     if( ( k = *(++B2) ) == InINF )
      break;
     if( ( h = *(++B1) ) == InINF )
      break;
     }

 }  // end( VectScale( g1 , B1 , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 *= x * g2

template<class T1, class T2, class T3>
inline void VectScale( register T1 *g1 , register const T2 *g2 , 
		       register const T3 x , register Index n )
{
 for( ; n-- ; )
  *(g1++) *= x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 *= x * g2{B}

template<class T1, class T2, class T3>
inline void VectScale( register T1 *g1 , register const T2 *g2 , 
		       register const T3 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g1++) *= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} *= x * g2{B}, B = intersection of B1 and B2, all other unchanged

template<class T1, class T2, class T3>
inline void VectScale( register T1 *g1 , register cIndex_Set B1 ,
		       register const T3 x , register const T2 *g2 ,
		       register cIndex_Set B2 )
{
 register Index h = *B1;
 register Index k = *B2;
 if( ( h < InINF ) && ( k < InINF ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == InINF )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == InINF )
      break;
     g1++;
     }
    else
    {
     (*(g1++)) *= x * (*(g2++));
     if( ( k = *(++B2) ) == InINF )
      break;
     if( ( h = *(++B1) ) == InINF )
      break;
     }

 }  // end( VectScale( g1 , B1 , x , g2 , B2 )

/*--------------------------------------------------------------------------*/
/// g1{B} *= g2, all other entries of g1 unchanged

template<class T1, class T2>
inline void VectScaleB( register T1 *g1 , register const T2 *g2 ,
			register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] *= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} *= x * g2, all other entries of g1 unchanged

template<class T1, class T2, class T3>
inline void VectScaleB( register T1 *g1 , register const T2 *g2 ,
			register const T3 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] *= x * (*(g2++));
 }

/*--------------------------------------------------------------------------*/
/// g1{B} *= g2{B}, all other entries unchanged

template<class T1, class T2>
inline void VectScaleBB( register T1 *g1 , register const T2 *g2 ,
			 register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] *= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} *= x * g2{B}, all other entries unchanged

template<class T1, class T2, class T3>
inline void VectScaleBB( register T1 *g1 , register const T2 *g2 ,
			 register const T3 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] *= x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// g := g / x, x a scalar
/** Divides each entry of g by the scalar x, in place. Useful e.g. for
 * *integer* types, for which VectScale( g , 1 / x , n ) would not work. */

template<class T>
inline void VectIScale( register T *g , register const T x ,
			register Index n )
{
 for( ; n-- ; )
  *(g++) /= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g{B} := g{B} / x, x a scalar

template<class T>
inline void VectIScale( register T *g , register const T x ,
		       register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g[ h ] /= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 /= g2

template<class T1, class T2>
inline void VectIScale( register T1 *g1 , register const T2 *g2 ,
		        register Index n )
{
 for( ; n-- ; )
  *(g1++) /= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 /= g2{B}

template<class T1, class T2>
inline void VectIScale( register T1 *g1 , register const T2 *g2 , 
		        register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g1++) /= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} /= g2{B}, B = intersection of B1 and B2, all other unchanged

template<class T1, class T2>
inline void VectIScale( register T1 *g1 , register cIndex_Set B1 ,
			register const T2 *g2 , register cIndex_Set B2 )
{
 register Index h = *B1;
 register Index k = *B2;
 if( ( h < InINF ) && ( k < InINF ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == InINF )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == InINF )
      break;
     g1++;
     }
    else
    {
     (*(g1++)) /= (*(g2++));
     if( ( k = *(++B2) ) == InINF )
      break;
     if( ( h = *(++B1) ) == InINF )
      break;
     }

 }  // end( VectIScale( g1 , B1 , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 /= x * g2

template<class T1, class T2, class T3>
inline void VectIScale( register T1 *g1 , register const T2 *g2 , 
		        register const T3 x , register Index n )
{
 for( ; n-- ; )
  *(g1++) /= ( x * (*(g2++)) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1 /= x * g2{B}

template<class T1, class T2, class T3>
inline void VectIScale( register T1 *g1 , register const T2 *g2 , 
		        register const T3 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g1++) /= ( x * g2[ h ] );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} /= x * g2{B}, B = intersection of B1 and B2, all other unchanged

template<class T1, class T2, class T3>
inline void VectIScale( register T1 *g1 , register cIndex_Set B1 ,
			register const T3 x , register const T2 *g2 ,
			register cIndex_Set B2 )
{
 register Index h = *B1;
 register Index k = *B2;
 if( ( h < InINF ) && ( k < InINF ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == InINF )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == InINF )
      break;
     g1++;
     }
    else
    {
     (*(g1++)) /= ( x * (*(g2++)) );
     if( ( k = *(++B2) ) == InINF )
      break;
     if( ( h = *(++B1) ) == InINF )
      break;
     }

 }  // end( VectIScale( g1 , B1 , x , g2 , B2 )

/*--------------------------------------------------------------------------*/
/// g1{B} /= g2, all other entries of g1 unchanged

template<class T1, class T2>
inline void VectIScaleB( register T1 *g1 , register const T2 *g2 ,
			 register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] /= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} /= x * g2, all other entries of g1 unchanged

template<class T1, class T2, class T3>
inline void VectIScaleB( register T1 *g1 , register const T2 *g2 ,
			 register const T3 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] /= x * (*(g2++));
 }

/*--------------------------------------------------------------------------*/
/// g1{B} /= g2{B}, all other entries unchanged

template<class T1, class T2>
inline void VectIScaleBB( register T1 *g1 , register const T2 *g2 ,
			  register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] *= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g1{B} /= x * g2{B}, all other entries unchanged

template<class T1, class T2, class T3>
inline void VectIScaleBB( register T1 *g1 , register const T2 *g2 ,
			  register const T3 x , register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] *= x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// g := g1 + x

template<class T1, class T2>
inline void VectAdd( register T1 *g , register const T2 *g1 ,
                     register const T1 x , register Index n )
{
 for( ; n-- ; )
  *(g++) = *(g1++) + x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := g1 + g2

template<class T1, class T2, class T3>
inline void VectAdd( register T1 *g , register const T2 *g1 ,
                     register const T3 *g2 , register Index n )
{
 for( ; n-- ; )
  *(g++) = *(g1++) + *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := g1 + g2{B}

template<class T1, class T2, class T3>
inline void VectAdd( register T1 *g , register const T2 *g1 ,
		     register const T3 *g2 , register cIndex_Set B ) 
{ 
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g++) = *(g1++) + g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := g1 + x * g2

template<class T1, class T2, class T3>
inline void VectAdd( register T1 *g , register const T2 *g1 ,
                     register const T3 *g2 , register const T3 x ,
                     register Index n )
{
 for( ; n-- ; )
  *(g++) = *(g1++) + x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := x1 * g1 + x

template<class T1, class T2>
inline void VectAdd( register T1 *g ,
		     register const T2 *g1 , register const T2 x1 ,
                     register const T1 x , register Index n )
{
 for( ; n-- ; )
  *(g++) = x1 * (*(g1++)) + x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := x1 * g1 + x2 * g2

template<class T1, class T2, class T3>
inline void VectAdd( register T1 *g ,
		     register const T2 *g1 , register const T2 x1 ,
                     register const T3 *g2 , register const T3 x2 ,
                     register Index n )
{
 for( ; n-- ; )
  *(g++) = x1 * (*(g1++)) + x2 * (*(g2++));
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// g := g1 - g2

template<class T1, class T2, class T3>
inline void VectDiff( register T1 *g , register const T2 *g1 ,
                      register const T3 *g2 , register Index n )
{
 for( ; n-- ; )
  *(g++) = *(g1++) - *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := g1 - g2{B}

template<class T1, class T2, class T3>
inline void VectDiff( register T1 *g , register const T2 *g1 ,
                      register const T3 *g2 , register cIndex_Set B ) 
{ 
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g++) = *(g1++) - g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// g := g1 * x

template<class T1, class T2>
inline void VectMult( register T1 *g , register const T2 *g1 ,
                      register const T1 x , register Index n )
{
 for( ; n-- ; )
  *(g++) = *(g1++) * x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := g1 * g2

template<class T1, class T2, class T3>
inline void VectMult( register T1 *g , register const T2 *g1 ,
                      register const T3 *g2 , register Index n )
{
 for( ; n-- ; )
  *(g++) = *(g1++) * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g += g1 * g2

template<class T1, class T2, class T3>
inline void VectMultAndSum( register T1 *g , register const T2 *g1 ,
                      register const T3 *g2 , register Index n )
{
 for( ; n-- ; )
  *(g++) += *(g1++) * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := g1 * g2{B}

template<class T1, class T2, class T3>
inline void VectMult( register T1 *g , register const T2 *g1 ,
                      register const T3 *g2 , register cIndex_Set B ) 
{ 
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g++) = *(g1++) * g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := x * g1 * g2

template<class T1, class T2, class T3, class T4>
inline void VectMult( register T1 *g ,
		      register const T2 x , register const T3 *g1 ,
                      register const T4 *g2 , register Index n )
{
 for( ; n-- ; )
  *(g++) = x * (*(g1++)) * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := x * g1 * g2{B}

template<class T1, class T2, class T3, class T4>
inline void VectMult( register T1 *g ,
		      register const T2 x , register const T3 *g1 ,
                      register const T4 *g2 , register cIndex_Set B ) 
{ 
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g++) = x * (*(g1++)) * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// g := g1 / g2

template<class T1, class T2, class T3>
inline void VectDivide( register T1 *g , register const T2 *g1 ,
		        register const T3 *g2 , register Index n )
{
 for( ; n-- ; )
  *(g++) = *(g1++) / (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := g1 / g2{B}

template<class T1, class T2, class T3>
inline void VectDivide( register T1 *g , register const T2 *g1 ,
		        register const T3 *g2 , register cIndex_Set B ) 
{ 
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g++) = *(g1++) / g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := x * ( g1 / g2 )

template<class T1, class T2, class T3, class T4>
inline void VectDivide( register T1 *g ,
		        register const T2 x , register const T3 *g1 ,
		        register const T4 *g2 , register Index n )
{
 for( ; n-- ; )
  *(g++) = x * ( (*(g1++)) / (*(g2++)) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g := x * ( g1 / g2{B} )

template<class T1, class T2, class T3, class T4>
inline void VectDivide( register T1 *g ,
		        register const T2 x , register const T3 *g1 ,
		        register const T4 *g2 , register cIndex_Set B ) 
{ 
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  *(g++) = x * ( (*(g1++)) / g2[ h ] );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// swap of g1 and g2 (element-wise)

template<class T>
inline void VectXcg( register T *g1 , register T *g2 , register Index n )
{
 for( ; n-- ; g1++ , g2++ )
  Swap( *g1 , *g2 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// swap of g1 and g2{B}

template<class T>
inline void VectXcg( register T *g1 , register T *g2 ,
                     register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; g1++ )
  Swap( *g1 , g2[ h ] );
 }

/*--------------------------------------------------------------------------*/
/// swap of g1{B} and g2

template<class T>
inline void VectXcgB( register T *g1 , register T *g2 ,
                     register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; g2++ )
  Swap( g1[ h ] , *g2 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// swap of g1{B} and g2{B}, all other entries are unchanged

template<class T>
inline void VectXcgBB( register T *g1 , register T *g2 ,
		       register cIndex_Set B )
{
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  Swap( g1[ h ] , g2[ h ] );
 }

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Functions for array manipulation (merging, shifting ...). The types  --*/
/*-- involved in these operations must support ordering relations '=='    --*/
/*-- and '>', as well as the assigment '='.                               --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Array manipulation
 *  @{ */

/// merges two ordered, "Stp-terminated" vectors into a third one
/** Merges the two ordered and "Stp-terminated" vectors g2 and g3 (that is,
 * there must be an element in both g2 and g3 that is >= Stp, which is
 * considered to be the terminator) into g, which will therefore also be
 * ordered and Stp-terminated.
 *
 * @param g     the vector where the merged, ordered, Stp-terminated result
 *              is written
 *
 * @param g2    the first ordered, Stp-terminated vector to be merged
 *
 * @param g3    the second ordered, Stp-terminated vector to be merged
 *
 * @param Stp   the common terminator value of g2 and g3. */

template<class T>
inline void Merge( register T *g , register const T *g2 ,
                   register const T *g3 , register const T Stp )
{
 register T i = *g2;
 register T j = *g3;

 for(;;)
  if( i < j )
  {
   if( ( *(g++) = i ) >= Stp )
    break;

   i = *(++g2);
   }
  else
  {
   if( ( *(g++) = j ) >= Stp )
    break;

   if( i == j )
    i = *(++g2);

   j = *(++g3);
   }
 }  // end( Merge )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// g[ 0 ] = g[ 1 ], g[ 1 ] = g[ 2 ] ... g[ n - 1 ] = g[ n ]

template<class T>
inline void ShiftVect( register T *g , register Index n )
{
 for( register T *t = g ; n-- ; t = g )
  *t = *(++g);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g[ i ] = g[ i + k ], i = 0 ... n - 1

template<class T>
inline void ShiftVect( register T *g , register Index n , cIndex k )
{
 for( register T *t = g + k ; n-- ; )
  *(g++) = *(t++);
 }

/*--------------------------------------------------------------------------*/
/// as ShiftVect(), but g[ 0 ] is moved to g[ n ]

template<class T>
inline void RotateVect( register T *g , register Index n )
{
 const T tmp = *g;

 for( register T *t = g ; n-- ; t = g )
  *t = *(++g);

 *g = tmp;
 }

/*--------------------------------------------------------------------------*/
/// g[ n - 1 ] = g[ n - 2 ], g[ n - 2 ] = g[ n - 3 ], ..., g[ 1 ] = g[ 0 ]

template<class T>
inline void ShiftRVect( register T *g , register Index n )
{
 for( register T *t = ( g += n ) ; n-- ; t = g )
  *t = *(--g);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// g[ i ] = g[ i - k ], i = n - 1 ... 0

template<class T>
inline void ShiftRVect( register T *g , register Index n , cIndex k )
{
 g += n;
 for( register T *t = g + k ; n-- ; )
  *(--t) = *(--g);
 }

/*--------------------------------------------------------------------------*/
/// as ShiftRVect(), but g[ n ] is moved to g[ 0 ]

template<class T>
inline void RotateRVect( register T *g , register Index n )
{
 register T *t = ( g += n );
 const T tmp = *t;

 for( ; n-- ; t = g )
  *t = *(--g);

 *g = tmp;
 }

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Searching/ordering functions. The types involved in these operations --*/
/*-- must support ordering relations '==' and '>', as well as the         --*/
/*-- assigment '='.                                                       --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Searching
 *  @{ */

/// returns the smallest index i such that g[ i ] == x, or n if none
/** If g contains elements identical to x, then the smallest index among all
 * such elements is reported (a number in 0 .. n - 1), else n is reported. */

template<class T>
inline Index Match( register const T *g , const T x , cIndex n )
{
 register Index i = 0;
 for( ; ( i < n ) && ( *(g++) != x ) ; )
  i++;

 return( i );
 }

/*--------------------------------------------------------------------------*/
/// returns TRUE <=> the two n-vectors g1 and g2 are element-wise identical

template<class T>
inline BOOL EqualVect( register const T *g1 ,  register const T *g2 ,
                       register Index n )
{
 for( ; n-- ; )
  if( *(g1++) != *(g2++) )
   return( FALSE );

 return( TRUE );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// as EqualVect( g1 , g2 , n ), but with g2 given in "sparse" form
/** Returns TRUE if and only if the two n-vectors g1 and g2 are element-wise
 * identical, where g2 is given in "sparse" form, i.e., g2[ i ] is the
 * B[ i ]-th element of the (implicitly "dense") vector to be compared with
 * g1, and all the entries of g1 whose index is not in B are required to be
 * == 0 for the vectors to be considered equal. */

template<class T>
inline BOOL EqualVect( register const T *g1 ,  register const T *g2 ,
                       register cIndex n , register cIndex_Set B )
{
 register Index i = 0;
 for( register Index h ; ( h = *(B++) ) < InINF ; i++ )
 {
  for( ; i < h ; i++ )
   if( *(g1++) )
    return( FALSE );

  if( *(g1++) != *(g2++) )
   return( FALSE );
  }

 for( ; i < n ; i++ )
  if( *(g1++) )
   return( FALSE );

 return( TRUE );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// binary search for What in the ordered, "infinity-terminated" set Set
/** Performs a binary search on the ordered set of Stop elements (of type T
 * and without replications) contained in the array Set: Set must be
 * "infinity-terminated", i.e., Set[ Stop ] > What. Searches for the
 * element What, that may or may not be in Set.
 *
 * @param Set    the ordered, "infinity-terminated" array to be searched
 *
 * @param Stop   the number of (valid) elements in Set
 *
 * @param What   the element being searched for
 *
 * @return if What is found, its position is reported; otherwise, the
 *         position of the smallest element > What in Set is reported. If
 *         there are no elements > What in Set, then Stop is reported. */

template<class T>
inline Index BinSearch( register const T *Set , register Index Stop ,
                        register const T What )
{
 register Index i = Stop / 2;
 for( register Index Strt = 0 ; Strt != Stop ; )
 {
  if( Set[ i ] == What )
   break;
  else
   if( Set[ i ] > What )
    Stop = i;
   else
    Strt = i + 1;

  i = ( Strt + Stop ) / 2;
  }

 return( i );

 }  // end( BinSearch )

/*--------------------------------------------------------------------------*/
/// like BinSearch() above, but What *must be* in Set (=> Set is nonempty)

template<class T>
inline Index BinSearch1( register const T *Set , register Index Stop ,
                         register const T What )
{
 register Index Strt = 0;
 register Index i = Stop / 2;
 for( register Index h ; ( h = Set[ i ] ) != What ; )
 {
  if( h > What )
   Stop = i - 1;
  else
   Strt = i + 1;

  i = ( Strt + Stop ) / 2;
  }

 return( i );

 }  // end( BinSearch1 )

/*--------------------------------------------------------------------------*/
/// like BinSearch(), but What is required *not* to be in Set
/** Like BinSearch(), but What is required *not* to be in Set, so that the
 * position of the smallest element > What in Set is always returned (Stop
 * if there is none). */

template<class T>
inline Index BinSearch2( register const T *Set , register Index Stop ,
                         register const T What )
{
 register Index i = Stop / 2;
 for( register Index Strt = 0 ; Strt != Stop ; )
 {
  if( Set[ i ] > What )
   Stop = i;
  else
   Strt = i + 1;

  i = ( Strt + Stop ) / 2;
  }

 return( i );

 }  // end( BinSearch2 )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/// inserts x into the binary heap H of n elements
/** H is a binary heap of elements of type T, ordered in increasing sense
 * (i.e., the root of the heap, H[ 1 ], is the smallest element) and
 * containing n elements: the new element x is inserted in H, which grows
 * to n + 1 elements. H is assumed to have room for the new element, i.e.,
 * H[ 1 .. n + 1 ] must be valid storage. */

template<class T>
inline void HeapIns( register T *H , const T x , register Index n )
{
 for( H-- , n++ ;;)
 {
  register Index p = n / 2;
  if( ( ! p ) || ( H[ p ] <= x ) )
   break;

  H[ n ] = H[ p ];
  n = p;
  }

 H[ n ] = x;

 }  // end( HeapIns )

/*--------------------------------------------------------------------------*/

/// deletes and returns the root (smallest element) of the binary heap H
/** H is a binary heap of elements of type T, as in HeapIns() [see above]
 * but here 0-based, i.e., H[ 0 ] is the root; returns the smallest element
 * (the root), deleting it from H and re-establishing the heap property.
 *
 * @param H   the (0-based) binary heap, modified in place
 *
 * @param n   the number of elements of H *after* the deletion, i.e.,
 *            | H | - 1 where | H | is the size of H before the call
 *
 * @return the smallest element of H (the root) prior to the deletion. */

template<class T>
inline Index HeapDel( register T *H , cIndex n )
{
 const T h = *H;
 *H = H[ n ];

 for( register Index i = 0 ; i < n  ; )
 {
  register Index j = 2 * i;

  if( ++j >= n )
   break;

  register Index k = j;

  if( ++k < n )
   if( H[ k ] < H[ j ] )
    k = j;

  if( H[ i ] <= H[ j ] )
   break;

  Swap( H[ i ] , H[ j ] );
  i = j;
  }

 return( h );
 }

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace OPTtypes_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* OPTvect.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File OPTvect.h --------------------------------*/
/*--------------------------------------------------------------------------*/
