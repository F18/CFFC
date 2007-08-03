/**********************************************************************
 * EmbeddedBoundaries2DInput.h: Header file defining the 2D embedded  *
 *                              boundary input parameter class.       *
 **********************************************************************/

#ifndef _EMBEDDEDBOUNDARIES2D_INPUT_INCLUDED
#define _EMBEDDEDBOUNDARIES2D_INPUT_INCLUDED

// Include the interface header file.

#ifndef _INTERFACE2D_INCLUDED
#include "Interface2D.h"
#endif // _INTERFACE2D_INCLUDED

// Define the structures and classes.

#define	INPUT_PARAMETER_LENGTH_EMBEDDEDBOUNDARIES2D  128

/*!
 * Class:  EmbeddedBoundaries2D_Input_Parameters
 *
 * @brief 
 *
 */
class EmbeddedBoundaries2D_Input_Parameters{
private:
public:
  //@{ @name Embedded boundary input parameters
  //! Interface counter.
  int ci;
  //! Number of interface components.
  int Number_of_Components;
  //! Interface input component list.
  Interface2D_List Component_List;
  //! Interface type name.
  char Type[INPUT_PARAMETER_LENGTH_EMBEDDEDBOUNDARIES2D];
  //! Interface boundary condition type.
  char BC_Type[INPUT_PARAMETER_LENGTH_EMBEDDEDBOUNDARIES2D];
  //! Interface motion type.
  char Motion[INPUT_PARAMETER_LENGTH_EMBEDDEDBOUNDARIES2D];
  //! Number of specified interface points.
  int *ns;
  //! Solution content redistribution type.
  char Solution_Redistribution_Type[INPUT_PARAMETER_LENGTH_EMBEDDEDBOUNDARIES2D];
  //! Solution content redistribution type.
  int Redistribution_Type;
  //! Interface evolution frequency.
  int Evolution_Frequency;
  //@}

  //@{ @name Constructors and destructor.
  //! Default constructor.
  EmbeddedBoundaries2D_Input_Parameters(void) {
    ci = 0;
    Number_of_Components = 0;
    strcpy(Type,"None");
    strcpy(BC_Type,"Reflection");
    strcpy(Motion,"Off");
    ns = NULL;
    strcpy(Solution_Redistribution_Type,"Injection");
    Redistribution_Type = SOLUTION_REDISTRIBUTION_INJECTION;
    Evolution_Frequency = 1;
  }

  //! Default destructor.
  ~EmbeddedBoundaries2D_Input_Parameters(void) {
    if (ns != NULL) { delete []ns; ns = NULL; }
    Component_List.deallocate();
  }
  //@}

  //@{ @name Broadcast routines
  void Broadcast_Input_Parameters(void);

#ifdef _MPI_VERSION
  void Broadcast_Input_Parameters(MPI::Intracomm &Communicator,
				  const int Source_CPU);
#endif
  //@}

  //@{ @name Input-output operators
  friend ostream &operator << (ostream &out_file,
		               const EmbeddedBoundaries2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       EmbeddedBoundaries2D_Input_Parameters &IP);
  //@}

};

/**********************************************************************
 * EmbeddedBoundaries2D_Input_Parameters::Broadcast_Input_Parameters  *
 **********************************************************************/
inline void EmbeddedBoundaries2D_Input_Parameters::Broadcast_Input_Parameters(void) {
#ifdef _MPI_VERSION
  // Interface input parameters:
  MPI::COMM_WORLD.Bcast(&(Number_of_Components),
			1,
			MPI::INT,0);
  // Broadcast the interface component list if necessary.
  if (Number_of_Components > 0)
    Broadcast_Interface_List(Component_List);
  // Solution content redistribution type:
  MPI::COMM_WORLD.Bcast(Solution_Redistribution_Type,
			INPUT_PARAMETER_LENGTH_EMBEDDEDBOUNDARIES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(Redistribution_Type),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(Evolution_Frequency),
			1,
			MPI::INT,0);
#endif
}

#ifdef _MPI_VERSION
inline void EmbeddedBoundaries2D_Input_Parameters::Broadcast_Input_Parameters(MPI::Intracomm &Communicator,
									      const int Source_CPU) {
  // Interface input parameters:
  MPI::COMM_WORLD.Bcast(&(Number_of_Components),
			1,
			MPI::INT,Source_CPU);
  // Broadcast the interface component list if necessary.
  if (Number_of_Components > 0)
    Broadcast_Interface_List(Component_List);
  // Solution content redistribution type:
  MPI::COMM_WORLD.Bcast(Solution_Redistribution_Type,
			INPUT_PARAMETER_LENGTH_EMBEDDEDBOUNDARIES2D,
			MPI::CHAR,Source_CPU);
  MPI::COMM_WORLD.Bcast(&(Redistribution_Type),
			1,
			MPI::INT,Source_CPU);
  MPI::COMM_WORLD.Bcast(&(Evolution_Frequency),
			1,
			MPI::INT,Source_CPU);
}
#endif

/**********************************************************************
 * EmbeddedBoundaries2D_Input_Parameters -- Input-output operators.   *
 **********************************************************************/
inline ostream &operator << (ostream &out_file,
			     const EmbeddedBoundaries2D_Input_Parameters &IP) {
  if (IP.Component_List.Ni) {
    out_file << "\n  -> Embedded Boundaries:";
    for (int ni = 1; ni <= IP.Component_List.Ni; ni++) {
      out_file << "\n     -> Interface " << ni << ": ";
      if (IP.Component_List[ni].Type == INTERFACE_CIRCLE) {
	out_file << "Circular interface";
      } else if (IP.Component_List[ni].Type == INTERFACE_ELLIPSE) {
	out_file << "Elliptical interface";
      } else if (IP.Component_List[ni].Type == INTERFACE_SQUARE) {
	out_file << "Square interface";
      } else if (IP.Component_List[ni].Type == INTERFACE_RECTANGLE) {
	out_file << "Rectangle interface";
      } else if (IP.Component_List[ni].Type == INTERFACE_ROCKET_PROPELLANT_GRAIN) {
	out_file << "SRM propellant grain";
      } else if (IP.Component_List[ni].Type == INTERFACE_NACA0012_AEROFOIL) {
	out_file << "NACA0012 aerofoil";
      } else if (IP.Component_List[ni].Type == INTERFACE_NACA0015_AEROFOIL) {
	out_file << "NACA0015 aerofoil";
      } else if (IP.Component_List[ni].Type == INTERFACE_NASA_ROTOR_37) {
 	//out_file << "NASA rotor 37 blade at " << IP.Rotor_Percent_Span << "% span";
 	out_file << "NASA rotor 37 blade";
      } else if (IP.Component_List[ni].Type == INTERFACE_NASA_ROTOR_67) {
 	//out_file << "NASA rotor 67 blade at " << IP.Rotor_Percent_Span << "% span";
 	out_file << "NASA rotor 67 blade";
      } else if (IP.Component_List[ni].Type == INTERFACE_ZALESAK) {
	out_file << "Zalesak's disk";
      } else if (IP.Component_List[ni].Type == INTERFACE_RINGLEB) {
	out_file << "Ringleb's flow domain";
      } else if (IP.Component_List[ni].Type == INTERFACE_BUMP) {
	out_file << "Bump channel flow";
      } else if (IP.Component_List[ni].Type == INTERFACE_FLAT_PLATE) {
	out_file << "Flat plate";
      } else if (IP.Component_List[ni].Type == INTERFACE_USER_SPECIFIED) {
	out_file << "User specified interface";
      }
      if (IP.Component_List[ni].Motion == INTERFACE_MOTION_STATIONARY) {
	out_file << " with no motion (stationary)";
      } else if (IP.Component_List[ni].Motion == INTERFACE_MOTION_EXPAND ||
		 IP.Component_List[ni].Motion == INTERFACE_MOTION_CONSTANT) {
	out_file << " with constant radial motion, vr = " 
		 << IP.Component_List[ni].Speed.x;
      } else if (IP.Component_List[ni].Motion == INTERFACE_MOTION_UNIFORM ||
		 IP.Component_List[ni].Motion == INTERFACE_MOTION_TRANSLATE) {
	out_file << " with constant translational motion";
      } else if (IP.Component_List[ni].Motion == INTERFACE_MOTION_ROTATE) {
	if (IP.Component_List[ni].Type == INTERFACE_CIRCLE) {
	  out_file << " with a periodic x-direction translation";
	} else if (IP.Component_List[ni].Type == INTERFACE_NACA0012_AEROFOIL ||
		   IP.Component_List[ni].Type == INTERFACE_NACA0015_AEROFOIL) {
	  out_file << " with a periodic oscillation about c/4";
	} else {
	  out_file << " with a constant rotation";
	}
      } else if (IP.Component_List[ni].Motion == INTERFACE_MOTION_BURNING_SURFACE) {
	out_file << " with a speed dependent on the local propellant burning rate";
      } else if (IP.Component_List[ni].Motion == INTERFACE_MOTION_MOMENTUM_TRANSFER) {
	out_file << " with motion dependent on the momentum transfer from the fluid";
      } else if (IP.Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STATIONARY) {
	out_file << " with no motion (stationary) from the solution of level set method";
      } else if (IP.Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_EXPAND) {
	out_file << " with constant radial motion, vr = " 
		 << IP.Component_List[ni].Speed.x << " from the solution of the level set method";
      } else if (IP.Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_STRETCH) {
	out_file << " with a stretched motion from the solution of the level set method";
      } else if (IP.Component_List[ni].Motion == INTERFACE_MOTION_LEVELSET_BULKFLOW) {
	out_file << " with motion dependent on the local flowfield (advection) from the solution of the level set method";
      }
    }
    out_file << "\n     -> Solution content redistribution type: " 
	     << IP.Solution_Redistribution_Type;
    out_file << "\n     -> Interface evolution frequency: " << IP.Evolution_Frequency;
  }

  return out_file;
}

inline istream &operator >> (istream &in_file,
			     EmbeddedBoundaries2D_Input_Parameters &IP) {
  return in_file;
}

#endif // _EMBEDDEDBOUNDARIES2D_INPUT_INCLUDED
