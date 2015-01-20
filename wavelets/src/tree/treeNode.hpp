
#ifndef TREENODE_H
#define TREENODE_H

#include <sstream>
#include "point.hpp"

template <typename T>
class TreeNode {
    public:
        explicit TreeNode(int j, int k, T inf, T sup, unsigned int nChilds, bool isFillable);
        TreeNode(const TreeNode<T> &other);
        TreeNode<T> & operator= (const TreeNode<T> &other);
        virtual ~TreeNode();

        int j() const;
        int level() const;

        int k() const;
        int offset() const;
        
        T inf() const;
        T sup() const;
        T center() const;

        Point<T> point() const;
        bool isFilled() const;
        bool isFillable() const;

        unsigned int nChilds() const;
        TreeNode** childs();
        TreeNode<T>* operator[](unsigned int childId) const;
        TreeNode<T>* operator()(T position) const;
        
        TreeNode<T>*& operator[](unsigned int childId);

        Point<T> setPoint(const T &point);
        
        virtual std::string toString() const;

        virtual TreeNode<T>* clone() const = 0;
        virtual void insert(const Point<T> &pt, unsigned int j = 0u) = 0;
        virtual int computeOffset(const Point<T> &pt, unsigned int j) = 0;

    protected:
        unsigned int _level;
        int _offset;

        T _inf, _sup, _center;
        
        unsigned int _nChilds;
        TreeNode **_childs;
    
        bool _isFillable, _isFilled;
        Point<T> _point;
};


template <typename T>
std::ostream & operator<< (std::ostream &os, const TreeNode<T> &node) {
    os << node.toString() << std::endl;
    return os;
}
        
template <typename T>
TreeNode<T>::TreeNode(int j, int k, T inf, T sup, unsigned int nChilds, bool isFillable) :
    _level(j), _offset(k),
    _inf(inf), _sup(sup), _center((_inf+_sup)/T(2)),
    _nChilds(nChilds), _childs(nullptr),
    _isFillable(isFillable), _isFilled(false), _point()
{
    assert(inf <= sup);

    _childs = new TreeNode<T>*[_nChilds];
    for (unsigned int i = 0; i < _nChilds; i++) {
        _childs[i] = nullptr;
    }
}

template <typename T>
TreeNode<T>::TreeNode(const TreeNode<T> &other) :
    _level(other.level()), _offset(other.offset()),
    _inf(other.inf()), _sup(other.sup()), _center((_inf+_sup)/T(2)),
    _nChilds(other.nChilds()), _childs(nullptr),
    _isFillable(other.isFillable()), _isFilled(other.isFilled()), _point(other.point())
{
    _childs = new TreeNode<T>*[_nChilds];
    for (unsigned int i = 0; i < this->nChilds(); i++) {
        _childs[i] = other[i]->clone();
    }
}
        
template <typename T>
TreeNode<T>& TreeNode<T>::operator= (const TreeNode<T> &other) {
    _level = other.level();
    _offset = other.offset();
    _inf = other.inf();
    _sup = other.sup();
    _center = (inf+sup)/T(2);
    _nChilds = other.nChilds();
    _isFillable = other.isFillable();
    _isFilled = other.isFilled();
    _point = other.point();

    delete [] _childs;
    _childs = nullptr;

    _childs = new TreeNode*[_nChilds];
    for (unsigned int i = 0; i < this->nChilds(); i++) {
        _childs[i] = other[i]->clone();
    }

    return *this;
}

template <typename T>
TreeNode<T>::~TreeNode() {
    for (unsigned int i = 0; i < _nChilds; i++) {
        delete _childs[i];
    }
    delete [] _childs;
}

template <typename T>
 int TreeNode<T>::j() const {
	 return _level;
}
template <typename T>
 int TreeNode<T>::level() const {
	 return _level;
}
template <typename T>
 int TreeNode<T>::k() const {
	 return _offset;
}
template <typename T>
 int TreeNode<T>::offset() const {
	 return _offset;
}
template <typename T>
 T TreeNode<T>::inf() const {
	 return _inf;
}
template <typename T>
 T TreeNode<T>::sup() const {
	 return _sup;
}
template <typename T>
 T TreeNode<T>::center() const {
	 return _center;
}
template <typename T>
 Point<T> TreeNode<T>::point() const {
	 return _point;
}
template <typename T>
 bool TreeNode<T>::isFilled() const {
	 return _isFilled;
}
template <typename T>
 bool TreeNode<T>::isFillable() const {
	 return _isFillable;
}
template <typename T>
unsigned int TreeNode<T>::nChilds() const {
	 return _nChilds;
}
template <typename T>
TreeNode<T>** TreeNode<T>::childs() {
    return _childs;
}
        
template <typename T>
TreeNode<T>* TreeNode<T>::operator[](unsigned int childId) const {
    assert(childId < _nChilds);
    return this->_childs[childId];
}

template <typename T>
TreeNode<T>*& TreeNode<T>::operator[](unsigned int childId) {
    assert(childId < _nChilds);
    return this->_childs[childId];
}

template <typename T>
TreeNode<T>* TreeNode<T>::operator()(T position) const {
    assert(position < sup && position > inf);
    unsigned int id = static_cast<unsigned int>((_nChilds-1)*(position-inf)/(sup-inf));
    return _childs[id];
}

template <typename T>
Point<T> TreeNode<T>::setPoint(const T &point) {
    T buffer = _point;
    _point = point;
    return buffer;
}

template <typename T>
std::string TreeNode<T>::toString() const {
    std::stringstream ss;
    ss << "TreeNode (" << this->level() << "," << this->offset() << ")" << std::endl;
    ss << "\tInterval [" << this->inf() << "," << this->sup() << "]" << std::endl;
    ss << "\tCenter x = " << this->center() << std::endl;
    ss << "\tFillable = " << this->isFillable() << "\t Filled = " << this->isFilled() << "\t Point = " << this->point() << std::endl;
    ss << "\t" << "Childs count : " << this->nChilds() << std::endl;
    return ss.str();
}

#endif /* end of include guard: TREENODE_H */
