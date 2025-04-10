#include <plplot.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define NSIZE    101
int main(int argc, char** argv) { 

  printf("number of args %d", argc);  

  (void) plparseopts( &argc, argv, PL_PARSE_FULL );
  PLFLT x[NSIZE], y[NSIZE];
  PLFLT xmin = 0., xmax = 1., ymin = 0., ymax = 100.;
  int   i;
  for ( i = 0; i < NSIZE; i++ )
  {
    x[i] = (PLFLT) ( i ) / (PLFLT) ( NSIZE - 1 );
    y[i] = ymax * x[i] * x[i] * x[i];
  }


  // Initialize plplot
  plinit();

  // Create a labelled box to hold the plot.
  plenv( xmin, xmax, ymin, ymax, 0, 0 );
  pllab( "x", "y=100 x#u3#d", "Simple PLplot demo of a 2D line plot" );

  // Plot the data that was prepared above.
  plline( NSIZE, x, y );

  // Close PLplot library
  plend();


  return EXIT_SUCCESS;

}
