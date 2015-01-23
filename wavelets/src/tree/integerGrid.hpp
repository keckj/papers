
#ifndef INTEGERGRID_H
#define INTEGERGRID_H

#include "treeNode.hpp"
#include "binaryTreeNode.hpp"

template <typename T>
class IntegerGrid : public TreeNode<T> {

    public:
        IntegerGrid(const Interval<T> &interval);
        IntegerGrid(const IntegerGrid<T> &other);
        ~IntegerGrid();
        
        virtual TreeNode<T>* clone() const override;
        virtual void makeChilds() override;
        virtual void insert(const Point<T> &pt, unsigned int j = 0u) override;
        virtual int computeOffset(const Point<T> &pt, unsigned int j) const override;
        virtual unsigned int computeChildId(T position) const override;

        int infInteger() const;
        int supInteger() const;
};

template <typename T>
IntegerGrid<T>::IntegerGrid(const Interval<T> &interval) :
    TreeNode<T>(-1,-1, interval.dilateToInteger(), nullptr, interval.integerCount(), false)
{
    this->makeChilds();
}

template <typename T>
IntegerGrid<T>::IntegerGrid(const IntegerGrid<T> &other) :
    TreeNode<T>(other) {
}

template <typename T>
IntegerGrid<T>::~IntegerGrid() {}
        
template <typename T>
TreeNode<T>* IntegerGrid<T>::clone() const {
    return new IntegerGrid<T>(*this);
}

template <typename T>
void IntegerGrid<T>::makeChilds() {

    //allocate even nodes
    for (int i = 0; i < this->nChilds(); i++) {
        this->child(i) = new BinaryTreeNode<T>(0, this->infInteger()+i, Interval<T>(this->infInteger()+i-1, this->infInteger()+i+1), this);
        if(i % 2 == 1) {
            this->child(i)->makeChilds();
            std::cout << i << std::endl;
        }
    }

    //allocate centered odd nodes
    for (unsigned int i = 1; i < this->nChilds() - 1u; i++) {
        if(i % 2 == 0) {
            this->child(i)->child(0) = this->child(i-1)->child(1);
            this->child(i)->child(1) = this->child(i+1)->child(0);
        }
    }

    //special border node handle
    this->child(0)->child(0) = nullptr;
    this->child(0)->child(1) = this->child(1)->child(0);

    if(this->nChilds() % 2 == 1)
        delete this->child(this->nChilds() - 1)->child(1);

    this->child(this->nChilds() - 1)->child(1) = nullptr;

    if(this->nChilds() % 2 == 1)
        this->child(this->nChilds() - 1)->child(0) = this->child(this->nChilds()-2)->child(1);
}

template <typename T>
void IntegerGrid<T>::insert(const Point<T> &pt, unsigned int j) {
    this->child(computeChildId(pt.x))->insert(pt, j);
}

template <typename T>
int IntegerGrid<T>::computeOffset(const Point<T> &pt, unsigned int j) const {
    T x = pt.x;
    assert(this->interval().contains(x));
    int id = static_cast<int>(x);
    //std::cout << x << " offset " << id << std::endl;
    return id;
}

template <typename T>
unsigned int IntegerGrid<T>::computeChildId(T position) const {
    return static_cast<unsigned int>(2*(position - this->inf())); 
}
        
template <typename T>
int IntegerGrid<T>::infInteger() const {
    return this->interval().lowerInteger();
}

template <typename T>
int IntegerGrid<T>::supInteger() const {
    return this->interval().upperInteger();
}

#endif
