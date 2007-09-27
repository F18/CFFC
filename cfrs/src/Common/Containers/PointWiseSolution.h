#ifndef _POINTWISESOLUTION_INCLUDED
#define _POINTWISESOLUTION_INCLUDED

#include "../../../../src_2D/Utilities/Utilities.h"
#include "Grid/Grid2D/Cell2D.h"
#include "Grid/Grid3D/Cell3D.h"

template <class NodeType = Node2D, class NodeSolutionType = double>
class PointWiseSolution;

/************************************************
*     Friend Functions : PointWiseSolution      *
************************************************/
template <class NodeType, class NodeSolutionType>
bool operator==(const PointWiseSolution<NodeType,NodeSolutionType>& left,
		const PointWiseSolution<NodeType,NodeSolutionType>& right);

template <class NodeType, class NodeSolutionType>
bool operator!=(const PointWiseSolution<NodeType,NodeSolutionType>& left,
		const PointWiseSolution<NodeType,NodeSolutionType>& right);

template <class NodeType, class NodeSolutionType>
std::ostream& operator<< (std::ostream& os, const PointWiseSolution<NodeType,NodeSolutionType>& Obj);

/*******************************************************
 * TEMPLETIZED CLASS: PointWiseSolution                *
 ******************************************************/
template <class NodeType, class NodeSolutionType>
class PointWiseSolution{
 private:
  NodeType Node;
  NodeSolutionType NodeSolution;
 public:

  // Default constructor
  PointWiseSolution(): NodeSolution(0.0){};
  // Constructor
  PointWiseSolution(const double & Val): NodeSolution(Val){ };
  PointWiseSolution(const NodeType &Node, const NodeSolutionType Solution):
    Node(Node), NodeSolution(Solution) { };
  // Destructor
  ~PointWiseSolution(){};
  // Copy Constructor
  PointWiseSolution( const PointWiseSolution &rhs): Node(rhs.Node), NodeSolution(rhs.NodeSolution){ };
  // Use the default assignment operator

  NodeSolutionType& GetValue(void) { return NodeSolution; }
  NodeType& GetNode(void) { return Node; }
  void SetValue(const NodeSolutionType & NodeSolution_) { NodeSolution = NodeSolution_;};
  void SetNode(const NodeType & Node_){ Node = Node_;};
  void SetParam(const NodeType & Node_ , const NodeSolutionType & NodeSolution_ ){
    SetValue(NodeSolution_); SetNode(Node_);
  }

  /* Friend functions */
  friend bool operator== <NodeType,NodeSolutionType> (const PointWiseSolution<NodeType,NodeSolutionType>& left,
						      const PointWiseSolution<NodeType,NodeSolutionType>& right);
  friend bool operator!= <NodeType,NodeSolutionType> (const PointWiseSolution<NodeType,NodeSolutionType>& left,
						      const PointWiseSolution<NodeType,NodeSolutionType>& right);
  /* output operator */
  friend std::ostream& operator<< <NodeType,NodeSolutionType> (std::ostream& os,
							       const PointWiseSolution<NodeType,NodeSolutionType> & pt);
};

/* Friend Functions */

/* Operator == */
template<class NodeType, class NodeSolutionType> inline
bool operator==(const PointWiseSolution<NodeType,NodeSolutionType>& left,
		const PointWiseSolution<NodeType,NodeSolutionType>& right){
  return (left.Node==right.Node)&&(left.NodeSolution==right.NodeSolution);
}

/* Operator != */
template<class NodeType, class NodeSolutionType> inline
bool operator!=(const PointWiseSolution<NodeType,NodeSolutionType>& left,
		const PointWiseSolution<NodeType,NodeSolutionType>& right){
  return !(left == right);
}

/* Operator << */
template<class NodeType, class NodeSolutionType> inline
std::ostream& operator<< (std::ostream& os,
			  const PointWiseSolution<NodeType,NodeSolutionType> & pt){
  os << "Node=" << pt.Node << " -> Solution=" << pt.NodeSolution << std::endl;
  return os;
}

#endif /*_POINTWISESOLUTION_INCLUDED */
