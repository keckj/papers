
#ifndef BINARYTREENODE_H
#define BINARYTREENODE_H

#include "treeNode.hpp"

template <typename T>
class BinaryTreeNode : public TreeNode<T> {

    public:
        explicit BinaryTreeNode(int j, int k, T inf, T sup);
        BinaryTreeNode(const BinaryTreeNode<T> &other);
        BinaryTreeNode<T> & operator= (const BinaryTreeNode<T> &other);
        ~BinaryTreeNode();
        
        virtual TreeNode<T>* clone() const override;
};
       
template <typename T>
BinaryTreeNode<T>::BinaryTreeNode(int j, int k, T inf, T sup) :
    TreeNode<T>(j,k,inf,sup,2u,true) {
}

template <typename T>
BinaryTreeNode<T>::BinaryTreeNode(const BinaryTreeNode<T> &other) : TreeNode<T>(other) {
}

template <typename T>
BinaryTreeNode<T>& BinaryTreeNode<T>::operator= (const BinaryTreeNode<T> &other) {
}

template <typename T>
BinaryTreeNode<T>::~BinaryTreeNode() {
}
        
template <typename T>
TreeNode<T>* BinaryTreeNode<T>::clone() const {
    return new BinaryTreeNode(*this);
}

#endif /* end of include guard: BINARYTREENODE_H */
