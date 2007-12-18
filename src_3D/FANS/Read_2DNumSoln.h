/**********************************************************************
 *Flowfield data for bluff-body stabilised turbulent non-reacting flows
 Burner bluff-body diameter, D_b=50mm, (R_b=25mm), jet diameter=3.6mm
 **********************************************************************/

#ifndef _READ2DNUMSOLN_INCLUDED
#define _READ2DNUMSOLN_INCLUDED
/* Include required C++ libraries. */

#include "Array/Array.h"

#ifndef _MULTIBLOCK_LIST_INCLUDED
#include "../HexaBlock/HexaMultiBlock.h"
#endif //_MULTIBLOCK_LIST_INCLUDED


/*******************************************************************************
 *
 * Routine: skip_line
 *
 * Purpose
 * =======
 *
 *   Discards a line.  The file position is set to the next character after the
 *   newline.  EOF is returned if it is detected.
 *
 ******************************************************************************/

inline int skip_line(FILE* fs)

{

   int a;
   do {
      a = fgetc(fs);
      if ( a == EOF ) return EOF;
   }
   while ( a != '\n' );
   return 0;
   
}


/*******************************************************************************
 *
 * Routine skip_field
 *
 * Purpose
 * =======
 *
 *   Discards all characters between file position and second set of whitespace.
 *   The file position is set to after the first whitespace character in the
 *   second set.  EOF is returned if detected.
 *
 ******************************************************************************/

inline int skip_field(FILE* fs)
   
{

   int a;
   do a = fgetc(fs);
   while ( std::isspace(a) );
   if ( a == EOF ) return EOF;
   do {
      a = fgetc(fs);
      if ( a == EOF ) return EOF;
   }
   while ( !isspace(a) );
   return 0;
   
}

class NumericalFlowField_2D {

  public:

   enum 
      {
         nzone = 43,
         ni = 12,
         nj = 12,
         nvar = 27
      };
   
   
   Array::Dynamic<double, 4> data;
   
   /* Creation constructors. */

   NumericalFlowField_2D(void): data(nzone, ni, nj, nvar){ }
   
   int read_numerical_solution(const char *const cffc_path);
   
   
};


inline int  NumericalFlowField_2D::read_numerical_solution(const char *const cffc_path){

   FILE * fp;

   char filename[256];
   
   strcpy(filename, cffc_path);
   strcat(filename, "/src_3D/InitialSolutionFields/BluffBodyBurner_hotflow_2D.dat");
   
   fp = fopen (filename,"r");
   
   if ( !fp ) {
      cout << "File not found.\n";
      return 1;
   }
   
   char str, zoneT;  
   int izone;
   
   for(int iblk = 0; iblk<nzone; iblk++){
      
      fscanf(fp, "%c", &str);
      if ( str == 'T' )
         for ( int n = 46 ; n-- ; ) skip_line(fp);  // Skip header lines
   
      // Debug
      fpos_t *fptr;
      fgetpos(fp, fptr);
      fgets(filename, 255, fp);
      
      std::cout << "Expecting zone:" << filename << std::endl;
      fsetpos(fp, fptr);
      // End debug
      fscanf(fp, "ZONE T =  \"Block Number = %d\"", &izone);
      
      cout<<"\n izone = "<<izone<<endl;
      
      skip_line(fp);
      skip_line(fp);
      skip_line(fp);
      skip_line(fp);
   
      for(int id = 0; id<ni; id++) 
         for(int jd = 0; jd<ni; jd++)
            for(int iv = 0; iv<nvar; iv++){
            
               fscanf(fp, "%lf", &data(izone, id, jd, iv));
               skip_line(fp);
            
            }


   }// finish reading all the data

   

   fclose (fp);
  


   return 0;
   
    
}


template <class HEXA_BLOCK>
int Read_2D_NumSoln(HEXA_BLOCK *Solution_Block, 
                    const int ICtypes, 
                    const char *const cffc_path) {
   
   if (Solution_Block[0].Flow_Type != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return (0);
   if (ICtypes == IC_TURBULENT_DIFFUSION_FLAME && (CFFC_Primary_MPI_Processor()) )

   
   return (0);
   
}




#endif /* _BLUFFBODYBURNER_INCLUDED  */


