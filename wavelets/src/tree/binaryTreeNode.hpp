
#ifndef BINARYTREENODE_H
#define BINARYTREENODE_H

#include "treeNode.hpp"

template <typename T>
class BinaryTreeNode : public TreeNode<T> {

    public:
        explicit BinaryTreeNode(int j, int k, Interval<T> interval, const TreeNode<T> *father);
        BinaryTreeNode(const BinaryTreeNode<T> &other);
        BinaryTreeNode<T> & operator= (const BinaryTreeNode<T> &other);
        ~BinaryTreeNode();

        bool canBeFather(const Point<T> &pt, unsigned int j, int k);
        
        virtual TreeNode<T>* clone() const override;
        virtual void insert(const Point<T> &pt, unsigned int j = 0u) override;
        virtual void makeChilds() override;
        virtual int computeOffset(const Point<T> &pt, unsigned int j) const override;
        virtual unsigned int computeChildId(T position) const override;
};
       
template <typename T>
BinaryTreeNode<T>::BinaryTreeNode(int j, int k, Interval<T> interval, const TreeNode<T>* father) :
    TreeNode<T>(j,k,interval,father,2u,true) {
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
    assert(this->_interval.contains(pt.x));

    //std::cout << "Tree node :  j:" << this->level() << "\tk:" << this->offset() << std::endl;
    //std::cout << "Point :  j:" << j << "\tk:" << k << "\tk/2**(j-l): " << (k/static_cast<int>(std::pow(2,j - this->level())))  << std::endl;

    return ( ((this->level() <= j) && (this->offset() == k))
        || ((this->level() < j) && 
             (k/static_cast<int>(std::pow(2,j - this->level()))) == this->offset()) );
}
        
template <typename T>
void BinaryTreeNode<T>::insert(const Point<T> &pt, unsigned int j) {

    int k = computeOffset(pt, j);
    assert(this->canBeFather(pt,j,k));
    
    //std::cout << "inserting in (" << j << "," << k <<  ") point " << pt.x << std::endl;
    //std::cout << "current node is (" << this->level() << "," << this->offset() << ") centered at " << this->center() << std::endl;
    
    T position = pt.x;
    TreeNode<T> *child;

    //tree descent
    if(this->level() < j) {
        this->makeChilds();
        child = this->getChild(computeChildId(position));
        child->insert(pt, j);
    }
    //insertion
    else if(this->level() == j) {
        this->makeChilds();

        //competition
        if(this->isFilled()) {
            T d = this->distanceFromCenter(this->position());
            T dnew = this->distanceFromCenter(position);

            //new point is better
            if(dnew < d) {
                Point<T> old = this->setPoint(pt);
                child = this->getChild(computeChildId(old.x));
                child->insert(old, j+1); // reinsert old point
            }
            // already affected point is better
            else {
                child = this->getChild(computeChildId(position));
                child->insert(pt, j+1);
            }
        }
        //affectation
        else {
            this->setPoint(pt);
        }
    }
    else {
        assert(false);
    }
}
        
template <typename T>
void BinaryTreeNode<T>::makeChilds() {
    
    if (this->_childs[0] == nullptr) {
        Interval<T> I1( this->inf(),    this->center() );
        Interval<T> I2( this->center(), this->sup()    );
        int J = this->level() + 1u;
        int K = this->offset();

        this->_childs[0] = new BinaryTreeNode<T>(J, 2*K-1, I1, this);
        this->_childs[1] = new BinaryTreeNode<T>(J, 2*K+1, I2, this);
    }
}
        
template <typename T>
unsigned int BinaryTreeNode<T>::computeChildId(T position) const {
    unsigned int childId;

    if(position < this->center())     
        childId = 0u;
    else
        childId = 1u;

    return childId;
}

template <typename T>
int BinaryTreeNode<T>::computeOffset(const Point<T> &pt, unsigned int j) const {

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
