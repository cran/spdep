/* Copyright 2001 by Nicholas Lewin-Koh. ****
* NOTE
* The subroutines TwoCirclesxx and SubVec are adapted for R and Double
* precision coordinates by Nicholas Lewin-Koh, from Computational
* Geometry in C, Joseph O.Rourke, Cambridge University Press
* (1998). Copyright for those subroutines remains his.
****/

#include <R.h>
#include <Rmath.h>


#define	X 0
#define	Y 1
#define DIM     2               /* Dimension of points */
typedef int     tPointi[DIM];   /* type integer point */
typedef double  tPointd[DIM];   /* type double point */

int     TwoCircles( double x1, double y1, double r1, 
		    double x2, double y2, double r2, tPointd p);
int	TwoCircles0a( double r1, tPointd c2, double r2, tPointd p );
int	TwoCircles0b( double r1, tPointd c2, double r2, tPointd p );
void	TwoCircles00( double r1, double a2, double r2, tPointd p );
void	SubVec( tPointd a, tPointd b, tPointd c );


static double distance(double x1, double y1, double x2, double y2){

  return(pythag(x1-x2,y1-y2));

}
static double  Length2( tPointd v )
{
   int   i;
   double   ss;

   ss = 0;
   for ( i=0; i < DIM; i++ )
      ss += R_pow_di(v[i],2);
   return ss;
}

/*---------------------------------------------------------------------
TwoCircles finds an intersection point between two circles.
General routine: no assumptions. Returns # of intersections; point in p.
---------------------------------------------------------------------*/
int     TwoCircles( double x1, double y1, double r1, 
		    double x2, double y2, double r2, tPointd p)
{
   tPointd c, q, c1, c2;
   int nsoln = -1;

   c1[X]=x1; 
   c1[Y]=y1;
   c2[X]=x2; 
   c2[Y]=y2;

   /* Translate so that c1={0,0}. */
   SubVec( c2, c1, c );
   nsoln = TwoCircles0a( r1, c, r2, q );
   /* Translate back. */
   p[X] = q[X] + c1[X];
   p[Y] = q[Y] + c1[Y];
   return nsoln;
}

/*---------------------------------------------------------------------
TwoCircles0a assumes that the first circle is centered on the origin.
Returns # of intersections: 0, 1, 2, 3 (inf); point in p.
---------------------------------------------------------------------*/
int     TwoCircles0a( double r1, tPointd c2, double r2, tPointd p )
{
   double dc2;              /* dist to center 2 squared */
   double rplus2, rminus2;  /* (r1 +/- r2)^2 */
   double f;                /* fraction along c2 for nsoln=1 */

   /* Handle special cases. */
   dc2 = Length2( c2 );
   rplus2  = R_pow_di((r1 + r2),2);
   rminus2 = R_pow_di((r1 - r2),2); 

   /* No solution if c2 out of reach + or -. */
   if ( ( dc2 > rplus2 ) || ( dc2 < rminus2 ) )
      return   0;

   /* One solution if c2 just reached. */
   /* Then solution is r1-of-the-way (f) to c2. */
   if ( dc2 == rplus2 ) {
      f = r1 / (r1 + r2);
      p[X] = f * c2[X];   p[Y] = f * c2[Y];
      return 1;
   }
   if ( dc2 == rminus2 ) {
      if ( rminus2 == 0 ) {   /* Circles coincide. */
         p[X] = r1;    p[Y] = 0;
         return 3;
      }
      f = r1 / (r1 - r2);
      p[X] = f * c2[X];   p[Y] = f * c2[Y];
      return 1;
   }

   /* Two intersections. */
   return TwoCircles0b( r1, c2, r2, p );
}

/*---------------------------------------------------------------------
TwoCircles0b also assumes that the 1st circle is origin-centered.
---------------------------------------------------------------------*/
int     TwoCircles0b( double r1, tPointd c2, double r2, tPointd p )
{
   double a2;          /* center of 2nd circle when rotated to x-axis */
   tPointd q;          /* one solution when c2 on x-axis */
   double cost, sint;  /* sine and cosine of angle of c2 */

   /* Rotate c2 to a2 on x-axis. */
   /*a2 = sqrt( Length2( c2 ) )*/ 
   a2=pythag(c2[X],c2[Y]);
   cost = c2[X] / a2;
   sint = c2[Y] / a2;

   TwoCircles00( r1, a2, r2, q );

   /* Rotate back */
   p[X] =  cost * q[X] + -sint * q[Y];
   p[Y] =  sint * q[X] +  cost * q[Y];
   
   return 2;
}

/*---------------------------------------------------------------------
TwoCircles00 assumes circle centers are (0,0) and (a2,0).
---------------------------------------------------------------------*/
void     TwoCircles00( double r1, double a2, double r2, tPointd p )
{
   double r1sq, r2sq;
   r1sq = R_pow_di(r1,2);
   r2sq = R_pow_di(r2,2);

   /* Return only positive-y soln in p. */
   p[X] = ( a2 + ( r1sq - r2sq ) / a2 ) * 0.5;
   p[Y] = pythag(r1,p[X]);
   /*sqrt( r1sq - p[X]*p[X] );*/
   /*Rprintf("%%TwoCircles00: p=(%f,%f)\n", p[X], p[Y]);*/
   return;
}

void   SubVec( tPointd a, tPointd b, tPointd c )
{
   int   i;

   for ( i=0; i < DIM; i++ ) {
      c[i] = a[i] - b[i];
   }
   return;
}
/* 
   This is a brute force intial implementation that will need a lot of
   refining. This computes the sphere of influence graph of a 2-d
   point set, using a neighbor list from a delaunay triangulation as
   input. The program makes 2 passes through the data.
 
*/ 

void compute_soi(int *no_nodes, int *g1, int *g2, int *noedges,
                  int *noneigh, int *neigh, int *nearneigh, 
                 double *rad, double *nodes_xd, double *nodes_yd)
 {
   int i, j, h=0, m=0;
   double temp;
   tPointd xtra;

   /* 
      First calculate the nearest neighbors of each point, and save
      the distance as the cricle radii for the next step.
   */ 
   for(i=0;i<*no_nodes;i++){
     for(j=0;j<noneigh[i];j++){
       if(rad[i]==0.0){
	 rad[i]=distance(nodes_xd[i], nodes_yd[i], nodes_xd[neigh[h]-1],
			 nodes_yd[neigh[h]-1]);
	 nearneigh[i]=neigh[h];
       }
       else{
	 temp=distance(nodes_xd[i], nodes_yd[i], 
		       nodes_xd[neigh[h]-1], nodes_yd[neigh[h]-1]);
	 if(temp < rad[i]){
	   rad[i]=temp;
	   nearneigh[i]=neigh[h];
	 }
       }
       h++;
     }
   }

   /* 
      Now for each point check for > 1 intersection of the points influence
      circle with all its Delaunay neighbors. If the above condition is met
      add an edge connecting the two points.
   */ 
   h=0;
   m=0;
   for(i=0;i<*no_nodes;i++){
     for(j=0;j<noneigh[i];j++){
       if(TwoCircles(nodes_xd[i], nodes_yd[i], rad[i], 
		     nodes_xd[neigh[h]-1], nodes_yd[neigh[h]-1], 
		     rad[neigh[h]-1],xtra) >= 2){
	 g1[m]=i+1;
	 g2[m]=neigh[h];
	 m++;
       }
       h++;
     }
   }
   *noedges=m;
   return;
 }

