#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

using namespace std;

/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 double first, second, diff;
 diff = -100000;

 string f1 = argv[ 1 ];
 string f2 = argv[ 2 ];

 ifstream iFile1( f1 );
 ifstream iFile2( f2 );
 int i=0;
 
 while (!iFile1.eof()) {
   i++;
   iFile1 >> first;
   iFile2 >> second;
   //cout << first-second << "\n";
   if( diff < abs( first - second ) && second < 10e6 )
     diff = abs( first - second );
   if( abs( first - second ) / first > 0.0001 ){
    std::cout << first << " " << second << "\n";
    //break;
   }
   
  }
 cout << diff << "\n";
  
 iFile1.close();
 iFile2.close();
  
 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
