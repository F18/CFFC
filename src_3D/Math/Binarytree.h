/*
 *  Binarytree.h
 *
 *  Created by Willem Deconinck on 11/02/08.
 *
 */
#ifndef _BINARY_TREE_INCLUDED
#define _BINARY_TREE_INCLUDED
#include <iostream>
using namespace std;

template<class itemtype>
class treeNode {
public:
    itemtype item;
    treeNode *left;
    treeNode *right;
    treeNode(itemtype &newItem) {
        item = newItem;
        left = NULL;
        right = NULL;
    }
};

template<class itemtype>
class BinaryTree {
protected:
    treeNode<itemtype> *root;     

public:
    
    // Default Constructor
    BinaryTree() : root(NULL) {  } 
    
    // Copy Constructor
    BinaryTree(const BinaryTree<itemtype>& otherTree); 
    
    // Destructor
    ~BinaryTree() { destroy(root); }
    
    //Overload the assignment operator
    const BinaryTree<itemtype>& operator= (const BinaryTree<itemtype>&); 
    
    
    void InsertNode(itemtype &newItem);
    void List( );
    int countNodes( );
    int treeHeight();
    itemtype *asArray();
    bool isEmpty();
    itemtype max( );
    itemtype min( );
    bool Contains(itemtype &item ) ;


    
    
private:
    int arrayIndex;
    void inorderTraversal(treeNode<itemtype> *node, itemtype *itemArray);
    void addToArray(itemtype *itemArray,itemtype item);
    void InsertNode(treeNode<itemtype> *&node, itemtype &newItem);
    bool Contains(treeNode<itemtype> *node, itemtype &item ) ;
    void List(treeNode<itemtype> *node);
    int countNodes(treeNode<itemtype> *node);
    void copyTree (treeNode<itemtype>* &copiedTreeRoot, treeNode<itemtype>* otherTreeRoot);
    void destroy(treeNode<itemtype>* &node);
    int height(treeNode<itemtype> *node);
    int max(int x, int y);
    itemtype max(treeNode<itemtype> *node);
    itemtype min(treeNode<itemtype> *node);

};        


template<class itemtype>
itemtype BinaryTree<itemtype>::max( ) {
    return max(root);
}


template<class itemtype>
itemtype BinaryTree<itemtype>::max(treeNode<itemtype> *node) {
    if(node->right == NULL)
        return node->item;
    else 
        return max(node->right); 
}


template<class itemtype>
itemtype BinaryTree<itemtype>::min( ) {
    return min(root);
}


template<class itemtype>
itemtype BinaryTree<itemtype>::min(treeNode<itemtype> *node) {
    if(node->left == NULL)
        return node->item;
    else 
        return min(node->left); 
}

template<class itemtype>
int BinaryTree<itemtype>::height(treeNode<itemtype> *node) {
	if(node == NULL)
		return 0;
	else
		return 1 + max(height(node->left), height(node->right));
}

template<class itemtype>
int BinaryTree<itemtype>::max(int x, int y) {
	if(x >= y)
		return x;
	else
		return y;
}

template<class itemtype>
bool BinaryTree<itemtype>::isEmpty()
{
	return (root == NULL);
}

template<class itemtype>
int BinaryTree<itemtype>::treeHeight()
{
	return height(root);
}

template<class itemtype>
void BinaryTree<itemtype>::InsertNode(itemtype &newItem) {
    InsertNode(root, newItem);
}
    
template<class itemtype>
void BinaryTree<itemtype>::InsertNode(treeNode<itemtype> *&node, itemtype &newItem) {
    if ( node == NULL ) {
        node = new treeNode<itemtype>( newItem );
        return;
    }
    else if ( newItem < root->item ) {
        InsertNode( node->left, newItem );
    }
    else {
        InsertNode( node->right, newItem );
    }
}
   

template<class itemtype>
bool BinaryTree<itemtype>::Contains(itemtype &item ) {
    if ( root == NULL ) {
        return false;
    }
    else if ( item == root->item ) {
        return true;
    }
    else if ( item < root->item ) {
        return Contains( root->left, item );
    }
    else {
        return Contains( root->right, item );
    }
}
        
template<class itemtype>
bool BinaryTree<itemtype>::Contains(treeNode<itemtype> *node, itemtype &item ) {
    if ( node == NULL ) {
        return false;
    }
    else if ( item == node->item ) {
        return true;
    }
    else if ( item < node->item ) {
        return Contains( node->left, item );
    }
    else {
        return Contains( node->right, item );
    }
}
    
template<class itemtype>
void BinaryTree<itemtype>::List( ) {
    List(root);
}
    
template<class itemtype>
void BinaryTree<itemtype>::List(treeNode<itemtype> *node) {
    // inorderTraversal
    if ( node != NULL ) {
        List(node->left);        
        cout << "  " << node->item << endl; 
        List(node->right); 
    }
}
    
template<class itemtype>
int BinaryTree<itemtype>::countNodes( ) {
    // Count the nodes in the binary tree to which node 
    // points.  Return the answer.
    if ( root == NULL ) {
        return 0;
    }
    else {
        // postorderTraversal
        int leftCount = countNodes( root->left );
        int rightCount = countNodes( root->right );
        return  1 + leftCount + rightCount;  
    }
} // end countNodes()
    
template<class itemtype>
int BinaryTree<itemtype>::countNodes(treeNode<itemtype> *node) {
    // Count the nodes in the binary tree to which node 
    // points.  Return the answer.
    if ( node == NULL ) {
        return 0;
    }
    else {
        // Add up the root node and the nodes in its two subtrees.
        int leftCount = countNodes( node->left );
        int rightCount = countNodes( node->right );
        return  1 + leftCount + rightCount;  
    }
} // end countNodes()


template<class itemtype>
itemtype * BinaryTree<itemtype>::asArray(){
    arrayIndex = 0;
    itemtype *itemlist = new itemtype[countNodes()];
    inorderTraversal(root, itemlist);
    return itemlist;
}


template<class itemtype>
void BinaryTree<itemtype>::inorderTraversal(treeNode<itemtype> *node, itemtype *itemArray){
    if ( node != NULL ) {
        inorderTraversal(node->left  , itemArray); // go to left most node         
        addToArray(itemArray, node->item); // start writing
        inorderTraversal(node->right , itemArray); // finish on the right
    }
}

template<class itemtype>
void BinaryTree<itemtype>::addToArray(itemtype *itemArray,itemtype item){
    itemArray[arrayIndex] = item;
    arrayIndex++;
}


template <class itemtype>
void  BinaryTree<itemtype>::destroy(treeNode<itemtype>* &node) {
	if(node != NULL) {
		destroy(node->left);
		destroy(node->right);
		delete node;
		node = NULL;
	}
}

//copy constructor
template <class itemtype>
BinaryTree<itemtype>::BinaryTree (const BinaryTree<itemtype>& otherTree) {
	if(otherTree.root == NULL) //otherTree is empty
		root = NULL;
	else
		copyTree(root, otherTree.root);
}


template <class itemtype>
void  BinaryTree<itemtype>::copyTree (treeNode<itemtype>* &copiedTreeRoot, treeNode<itemtype>* otherTreeRoot)
{
	if(otherTreeRoot == NULL)
		copiedTreeRoot = NULL;
	else
	{
		copiedTreeRoot = new treeNode<itemtype>;
		copiedTreeRoot->item = otherTreeRoot->item;
		copyTree(copiedTreeRoot->left, otherTreeRoot->left);
		copyTree(copiedTreeRoot->right, otherTreeRoot->right);
	}
} //end copyTree

template<class itemtype>
const BinaryTree<itemtype>& BinaryTree<itemtype>:: operator=(const BinaryTree<itemtype>& otherTree) { 
    
	if(this != &otherTree) { //avoid self-copy
		if(root != NULL)  //if the binary tree is not empty, destroy the binary tree
			destroy(root);
        
		if(otherTree.root == NULL) //otherTree is empty
			root = NULL;
		else
			copyTree(root, otherTree.root);
	}
    
    return *this; 
}




#endif


