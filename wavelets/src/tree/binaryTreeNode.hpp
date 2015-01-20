
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

        bool canBeFather(const Point<T> &pt, unsigned int j, int k);
        
        virtual TreeNode<T>* clone() const override;
        virtual void insert(const Point<T> &pt, unsigned int j = 0u) override;
        virtual int computeOffset(const Point<T> &pt, unsigned int j) override;
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
        
template <typename T>
bool BinaryTreeNode<T>::canBeFather(const Point<T> &pt, unsigned int j, int k) {
    assert(pt.x >= this->inf() && pt.x <= this->sup());

    std::cout << "Tree node :  j:" << this->level() << "\tk:" << this->offset() << std::endl;
    std::cout << "Point :  j:" << j << "\tk:" << k << "\tk/2**(j-l): " << (k/static_cast<int>(std::pow(2,j - this->level())))  << std::endl;

    return ( ((this->level() <= j) && (this->offset() == k))
        || ((this->level() < j) && 
             (k/static_cast<int>(std::pow(2,j - this->level()))) == this->offset()) );
}
        
template <typename T>
void BinaryTreeNode<T>::insert(const Point<T> &pt, unsigned int j) {

    std::cout << "inserting in node centered in " << this->center() << std::endl;
    std::cout << "with point " << pt.x << std::endl;
    int k = computeOffset(pt, j);

    assert(this->canBeFather(pt,j,k));

    //tree descent
    if(this->level() < j) {
         
    }
}

template <typename T>
int BinaryTreeNode<T>::computeOffset(const Point<T> &pt, unsigned int j) {

    assert(j > 0);
    int k = floor(pt.x*std::pow(2, j-1));
    return 2*k+1;

    //T offset = T(1)/std::pow(T(2), T(j));
    //T x = pt.x;
    //std::cout << "x = " << x << std::endl;
    //std::cout << "offset = " << offset << std::endl;
    //std::cout << "k = " << k << std::endl;
    //std::cout << std::endl;
}

#endif /* end of include guard: BINARYTREENODE_H */
