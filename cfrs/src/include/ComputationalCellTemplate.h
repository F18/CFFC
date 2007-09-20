// ComputationalCellTemplate.h defines the base template for ComputationalCell class

#ifndef _COMPUTATIONALCELLTEMPLATE_INCLUDED
#define _COMPUTATIONALCELLTEMPLATE_INCLUDED

// Define the ComputationalCell class
template< SpaceType SpaceDimension = TwoD, class GeometryType = Cell2D_Cartesian, class SolutionType = double>
class ComputationalCell;

 /****************************************************************************
 * TEMPLATIZED CLASS: ComputationalCell                                     *
 *                                                                          *
 * VARIABLES:                                                               *
 *  CELL VARIABLES                                                          *
 *   Cell Geometry                                                          *
 *     geom      -- Stores the geometry of the cell.                        *
 *                  Possible types:                                         *
 *                  1D: Cell1D_Uniform, Cell1D_NonUniform                   *
 *                  2D: Cell2D, Cell2D_Cartesian, Cell2D_Quad               *
 *                                                                          *
 *   Cell Solution                                                          *
 *     U_cell    -- The average cell quantity. Type: SolutionType           *
 *     TD        -- TaylorDerivatives container.                            *
 *     GeomCoeff -- stores the geometric coefficients of the cell           *
 *                  These coeffs represent integrals of the form            *
 *                   I_nm= 1/A* Int((x-xCC)^n*(y-yCC)^m)                    *
 *     rings     -- Shows the number of layers around the reconstructed     *
 *                  cell that are used for the reconstruction in that cell  *
 *     RO        -- The desired Order of Reconstruction                     *
 *     FinalOrder -- The final Order of Reconstruction                      *
 *                                                                          *
 *  SUBDOMAIN VARIABLES                                                     *
 *      In each computational cell a subdomain grid is defined. At these    *
 *      points the solution is computed for assessing the accuracy.         *
 *                                                                          *
 *     NbSubgridPoints -- Number of subgrid points in the X, Y              *
 *                        and Z directions.                                 *
 *                                                                          *
 *    Solution variables in the subdomain                                   *
 *     SubGridSolution -- Solution in the subgrid domain                    *
 *                        Variable of type                                  *
 *                        SubGridMesh<Node,SpaceDimension,SolutionType>     *
 *                                                                          *
 *  ACCURACY                                                                *
 *     ErrorRec = average reconstruction error                              *
 *                ErrorRec = Integral(|F_numeric - F_analytic|)/ Area       *
 *                                                                          *
 ***************************************************************************/

/* *****************************************************************************************
********************************************************************************************
              END COMPUTATIONAL CELL BASE TEMPLATE CLASS
********************************************************************************************/

#endif
