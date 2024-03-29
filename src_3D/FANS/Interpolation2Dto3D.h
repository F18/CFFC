/* Interpolation2Dto3D.h: */  

#ifndef _INTERPOLATION2DTO3D_INCLUDED
#define _INTERPOLATION2DTO3D_INCLUDED

/* Include required CFFC header files. */

#include "Array/Array.h"

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

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


class FlowField_2D {

  public:

   int nzone;
   int ni;
   int nj;
   int nvar;
   

     
   
   Array::Dynamic<double, 4> data;
   
   /* Creation constructors. */
   FlowField_2D(void):nzone(43), ni(10), nj(10), nvar(30){};
   
   int read_numerical_solution_hotflow(const char *const cffc_path);
   int read_numerical_solution_coldflow(const char *const cffc_path);
   
   void Allocate(int i_ICs){
      if(i_ICs == IC_TURBULENT_DIFFUSION_FLAME){
         nzone = 43;
         ni = 12;
         nj = 12;
         nvar = 30;
         
         data.allocate(nzone, ni, nj, nvar);
      }
      if(i_ICs == IC_TURBULENT_COFLOW){
         
         nzone = 14;
         ni = 20;
         nj = 20;
         nvar = 30;
         data.allocate(nzone, ni, nj, nvar);
      }
   };
   

   void Deallocate(int i_ICs){
      if(i_ICs == IC_TURBULENT_DIFFUSION_FLAME){
         data.deallocate();
      }
      if(i_ICs == IC_TURBULENT_COFLOW){
         data.deallocate();
      }
      
   };
   
   
};


inline int  FlowField_2D::read_numerical_solution_hotflow(const char *const cffc_path){

      
   FILE * fp;
   
   char filename[256];
   
   strcpy(filename, cffc_path);
   strcat(filename, "/src_3D/FANS/BluffBodyBurner_hotflow_2D.dat");
   
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
      
      fscanf(fp, "ZONE T =  \"Block Number = %d \"", &izone);
      
      skip_line(fp);
      skip_line(fp);
      skip_line(fp);
      skip_line(fp);
   
     
      for(int jd = 0; jd<nj; jd++){
         for(int id = 0; id<ni; id++){
            for(int iv = 0; iv<nvar; iv++){
               
               fscanf(fp, "%lf", &data(izone, id, jd, iv));
              
               
            }
            skip_line(fp);
         }
      }
      
      
   }// finish reading all the data

   
   
   fclose (fp);
  


   return 0;
   
    
}

inline int  FlowField_2D::read_numerical_solution_coldflow(const char *const cffc_path){

   
   FILE * fp;

   char filename[256];
   
   strcpy(filename, cffc_path);
   strcat(filename, "/src_3D/FANS/BluffBodyBurner_coldflow_2D.dat");
   
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
         for ( int n = 37 ; n-- ; ) skip_line(fp);  // Skip header lines
      
      fscanf(fp, "ZONE T =  \"Block Number = %d \"", &izone);
      
      skip_line(fp);
      skip_line(fp);
      skip_line(fp);
      skip_line(fp);
   
     
      for(int jd = 0; jd<nj; jd++){
         for(int id = 0; id<ni; id++){
            for(int iv = 0; iv<nvar; iv++){
               
               fscanf(fp, "%lf", &data(izone, id, jd, iv));
              
               
            }
            skip_line(fp);
         }
      }
      
      
   }// finish reading all the data

   fclose (fp);
  
   return 0;
   
}

#endif // _INTERPOLATION2DTO3D_INCLUDED

