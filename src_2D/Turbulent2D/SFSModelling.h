/****************** SFSModelling.h **************************************
   This class defines some tools for eddy-viscosity type subfilter 
   scale models.
 ************************************************************************/
#ifndef _SFS_MODELLING_INCLUDED
#define _SFS_MODELLING_INCLUDED


#ifndef _TENSOR2D_INCLUDED
#include "../Math/Tensor2D.h"
#endif //_TENSOR2D_INCLUDED


/*!
 *
 * Class: SubfilterScaleModels
 *
 */
class SubfilterScaleModels{

  public:
    static double Smagorinsky_coef; //!< Smagorinsky coefficient.
    static double CI_Yoshizawa;     //!< Coefficient to calculate k based on Yoshizawa's model.
    static double CV_coef;          //!< Coefficient to estimate the eddy viscosity in the SFS k-equation.  
    static double CEPS_coef;        //!< Coefficient to calculate the dissipation rate of SFS k.


    // Two coefficients for Smagorinsky and Yoshizawa or two coefficients for k-equation.
    void set_SFSmodel_constants(const double &, const double &);
    void set_SFSmodel_constants_k(const double &, const double &);

    // sqrt(2*S*S)
    double abs_strain_rate(const Tensor2D &strain_rate) const;

    //! Kinematic turbulent viscosity using Smagorinsky model
    double eddy_viscosity_Smagorinsky(const Tensor2D &strain_rate,
				      const double &filter_width) const;

    //! Eddy viscosity based on k
    double eddy_viscosity_k(const double &k, 
			    const double &filter_width) const;

    //! k based on Yoshizawa's model 
    double sfs_k_Yoshizawa(const Tensor2D &strain_rate,
			   const double &filter_width) const;
       
};



#endif
