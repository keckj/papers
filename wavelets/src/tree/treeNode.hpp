
#ifndef TREENODE_H
#define TREENODE_H

#include <sstream>
#include <Eigen/Dense>
#include "point.hpp"
#include "gnuplot.hpp"
#include "waveletMapper.hpp"

template <typename T>
class TreeNode {
    public:
        explicit TreeNode(int j, int k, Interval<T> interval, const TreeNode<T> *father, unsigned int nChilds, bool isFillable);
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
        Interval<T> interval() const;

        Point<T> point() const;
        T position() const;
        T value() const;
        bool isFilled() const;
        bool isFillable() const;
        bool isValid() const;

        const TreeNode<T>* father() const;

        unsigned int nChilds() const;
        TreeNode<T>** childs();
        TreeNode<T>* getChild(unsigned int childId) const;
        TreeNode<T>* operator[](unsigned int childId) const;
        TreeNode<T>*& operator[](unsigned int childId);

        TreeNode<T>* operator()(T position) const;

        Point<T> setPoint(const Point<T> &point);
        static T distance(T p1, T p2);
        T distanceFromCenter(T position);

        unsigned int maxLevel() const;
        unsigned int countValidNodes() const;
        void getValidPositions(std::vector<T> &positions) const;
        void fillSystem(Eigen::MatrixXf &A, Eigen::VectorXf &b, const WaveletMapper<T> &waveletMapper) const;

        void plot(Gnuplot &gp) const;
        void plotValid(Gnuplot &gp) const;
        
        virtual std::string toString() const;

        virtual TreeNode<T>* clone() const = 0;
        virtual void makeChilds() = 0;
        virtual void insert(const Point<T> &pt, unsigned int j = 0u) = 0;
        virtual int computeOffset(const Point<T> &pt, unsigned int j) = 0;
        virtual unsigned int computeChildId(T position) = 0;
        
        void computePlacementCondition(unsigned int p, T rho);

    protected:
        unsigned int _level;
        int _offset;

        Interval<T> _interval;

        const TreeNode *_father;
        
        unsigned int _nChilds;
        TreeNode **_childs;
    
        bool _isFillable, _isFilled, _isValid;
        Point<T> _point;
        
        void setValidity(bool valid);

        void plotRec(std::vector<std::tuple<T,T>> &pts, unsigned int jmax) const;
        void plotValidRec(std::vector<std::tuple<T,T>> &pts, unsigned int jmax) const;
        
        unsigned int fillSystemRec(Eigen::MatrixXf &A, Eigen::VectorXf &b, 
                const WaveletMapper<T> &waveletMapper, 
                const std::vector<T> &validPositions, 
                unsigned int column) const;
        
        bool placementConditionUp(unsigned int p, T rho, const TreeNode<T> *checkedNode) const;
        bool placementConditionDown(unsigned int p, T rho, const TreeNode<T> *checkedNode);

        static bool checkFirstPlacementCondition(unsigned int p, T rho, const TreeNode<T> *node, const TreeNode<T> *checkedNode);
        static bool checkSecondPlacementCondition(T rho, const TreeNode<T> *checkedNode);
};


template <typename T>
std::ostream & operator<< (std::ostream &os, const TreeNode<T> &node) {
    os << node.toString() << std::endl;
    return os;
}
        
template <typename T>
TreeNode<T>::TreeNode(int j, int k, Interval<T> interval, const TreeNode<T> *father, unsigned int nChilds, bool isFillable) :
    _level(j), _offset(k),
    _interval(interval), 
    _father(father),
    _nChilds(nChilds), _childs(nullptr),
    _isFillable(isFillable), _isFilled(false), _isValid(false), _point()
{
    _childs = new TreeNode<T>*[_nChilds];
    for (unsigned int i = 0; i < _nChilds; i++) {
        _childs[i] = nullptr;
    }
}

template <typename T>
TreeNode<T>::TreeNode(const TreeNode<T> &other) :
    _level(other.level()), _offset(other.offset()),
    _interval(other.interval()),
    _father(other.father()),
    _nChilds(other.nChilds()), _childs(nullptr),
    _isFillable(other.isFillable()), _isFilled(other.isFilled()), _isValid(other.isValid()), _point(other.point())
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
    _interval = other.interval();
    _nChilds = other.nChilds();
    _father = const_cast<TreeNode<T>*>(other.father());
    _isFillable = other.isFillable();
    _isFilled = other.isFilled();
    _isValid = other.isValid();
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
	 return _interval.inf;
}
template <typename T>
 T TreeNode<T>::sup() const {
	 return _interval.sup;
}
template <typename T>
 T TreeNode<T>::center() const {
	 return _interval.center();
}
template <typename T>
Interval<T> TreeNode<T>::interval() const {
    return _interval;
}
template <typename T>
 Point<T> TreeNode<T>::point() const {
	 return _point;
}
template <typename T>
T TreeNode<T>::position() const {
    return _point.x;
}
template <typename T>
T TreeNode<T>::value() const {
    return _point.y;
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
 bool TreeNode<T>::isValid() const {
	 return _isValid;
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
TreeNode<T>* TreeNode<T>::getChild(unsigned int childId) const {
    return _childs[childId];
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
Point<T> TreeNode<T>::setPoint(const Point<T> &point) {
    Point<T> buffer = _point;
    _point = point;
    _isFilled = true;
    return buffer;
}

template <typename T>
T TreeNode<T>::distanceFromCenter(T position) {
    return sqrt((position - this->center())*(position - this->center()));
}
    
template <typename T>
T TreeNode<T>::distance(T p1, T p2) {
    return sqrt((p2-p1)*(p2-p1));
}
        
template <typename T>
unsigned int TreeNode<T>::maxLevel() const {

    unsigned int newlevel = 0u;
    TreeNode<T> *node;

    for (unsigned int i = 0; i < _nChilds; i++) {
        node = _childs[i];
        if(node != nullptr)
            newlevel = std::max(newlevel, 1u + node->maxLevel());
    }
    
    return newlevel;
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
        
template <typename T>
void TreeNode<T>::plot(Gnuplot &gp) const {
    unsigned int jmax = this->maxLevel();

    gp << "set xr [" << this->inf() << ":" << this->sup() <<  "]\n";
    gp << "set yr [0:" << jmax + 2 << "]\n";
    gp << "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5\n";

    std::vector<std::tuple<T,T>> pts;
    this->plotRec(pts, jmax);
    
    gp << "plot '-' with linespoints ls 1 title 'allocation'\n";
    gp.send1d(pts);
}

template <typename T>
void TreeNode<T>::plotRec(std::vector<std::tuple<T,T>> &pts, unsigned int jmax) const {
   
    std::tuple<T,T> c = std::make_tuple<T,T>(this->position(), jmax + 1u - this->level());

    //insert child points
    TreeNode<T> *node; 
    for (unsigned int i = 0; i < this->nChilds(); i++) {
        node = _childs[i];
        if(node != nullptr && node->isFilled()) {
            pts.push_back(c); 
            node->plotRec(pts, jmax);
        }
    }
    
    pts.push_back(c); 
}

template <typename T>
void TreeNode<T>::plotValid(Gnuplot &gp) const {
    unsigned int jmax = this->maxLevel();

    std::vector<std::tuple<T,T>> pts;
    this->plotValidRec(pts, jmax);
    
    gp << "plot '-' with point pt 7\n";
    gp.send1d(pts);
}

template <typename T>
void TreeNode<T>::plotValidRec(std::vector<std::tuple<T,T>> &pts, unsigned int jmax) const {
   
    std::tuple<T,T> c = std::make_tuple<T,T>(this->position(), jmax + 1u - this->level());

    //insert child points
    TreeNode<T> *node; 
    for (unsigned int i = 0; i < this->nChilds(); i++) {
        node = _childs[i];
        if(node != nullptr && node->isFilled()) {
            node->plotValidRec(pts, jmax);
        }
    }
   
    if(this->isValid())
        pts.push_back(c); 
}

template <typename T>
void TreeNode<T>::setValidity(bool valid) {
    this->_isValid = valid;
}

template <typename T>
void TreeNode<T>::computePlacementCondition(unsigned int p, T rho) {
    
    bool isValid = TreeNode<T>::checkSecondPlacementCondition(rho, this);
    isValid = isValid && this->placementConditionUp(p,rho,this);

    this->setValidity(isValid);

    TreeNode<T> *child;
    for (unsigned int i = 0; i < _nChilds; i++) {
        child = this->getChild(i);
        if(child != nullptr)
            child->computePlacementCondition(p,rho);
    }
}

template <typename T>
bool TreeNode<T>::placementConditionUp(unsigned int p, T rho, const TreeNode<T> *checkedNode) const {
    
    bool isValid = true;
    const TreeNode<T> *father = this->father();

    if(father == nullptr)
        return isValid;

    //check current node's father
    isValid = TreeNode<T>::checkFirstPlacementCondition(p, rho, checkedNode, father);
    if(!isValid)
        return false;

    //check current node's brothers and little brothers
    TreeNode<T> *brother;
    for (unsigned int i = 0; i < _nChilds; i++) {
        brother = father->getChild(i);
        //if the brother is not the current node and is filled
        if(brother != nullptr && brother->offset() != this->offset() && brother->isFilled()) { 
            isValid = brother->placementConditionDown(p, rho, checkedNode);
            if(!isValid)
                return false;
        }
    }

    //check current node's grand father
    isValid = father->placementConditionUp(p, rho, checkedNode);

    return isValid;
}

template <typename T>
bool TreeNode<T>::placementConditionDown(unsigned int p, T rho, const TreeNode<T> *checkedNode) {

    bool isValid;

    //current node's level is coarser or equal to the checked one !
    if(this->level() >= checkedNode->level())
        return true;

    //TODO speed up tree descent

    //check current node
    isValid = TreeNode<T>::checkFirstPlacementCondition(p, rho, checkedNode, this);
    if(!isValid)
        return false;
    
    //check current node's childs
    TreeNode<T> *child;
    for (unsigned int i = 0; i < _nChilds; i++) {
        child = this->getChild(i);
        //if the child is not the current node and is filled
        if(child != nullptr && child->offset() != this->offset() && child->isFilled()) { 
            isValid = child->placementConditionDown(p, rho, checkedNode);
            if(!isValid)
                return false;
        }
    }

    //all checked we can validate the current node
    return true;
}
        
template <typename T>
bool TreeNode<T>::checkFirstPlacementCondition(unsigned int p, T rho, const TreeNode<T> *node, const TreeNode<T> *checkedNode) {
    if(TreeNode<T>::distance(node->center(), checkedNode->center()) <= T(p)/std::pow(2,checkedNode->level())) {
        return (TreeNode<T>::distance(node->center(), node->position()) <= rho/std::pow(2,checkedNode->level())); 
    }
    else {
        return true;
    }
}

template <typename T>
bool TreeNode<T>::checkSecondPlacementCondition(T rho, const TreeNode<T> *checkedNode) {
    return (TreeNode<T>::distance(checkedNode->center(), checkedNode->position()) <= rho/std::pow(2,checkedNode->level()));
}
        
template <typename T>
const TreeNode<T>* TreeNode<T>::father() const {
    return _father;
}
        
template <typename T>
unsigned int TreeNode<T>::countValidNodes() const {

    unsigned int validNodes;

    if(this->isValid())
        validNodes = 1;
    else
        validNodes = 0;
    
    TreeNode<T> *child;
    for (unsigned int i = 0; i < _nChilds; i++) {
        child = _childs[i];
        if(child != nullptr)
            validNodes +=  child->countValidNodes();
    }
    
    return validNodes;
}
        
template <typename T>
void TreeNode<T>::getValidPositions(std::vector<T> &positions) const {

    if(this->isValid())
        positions.push_back(this->position());
    
    TreeNode<T> *child;
    for (unsigned int i = 0; i < _nChilds; i++) {
        child = _childs[i];
        if(child != nullptr)
            child->getValidPositions(positions);
    }
}

template <typename T>
void TreeNode<T>::fillSystem(Eigen::MatrixXf &A, Eigen::VectorXf &b, 
        const WaveletMapper<T> &waveletMapper) const {
    
    std::vector<T> validPositions;

    this->getValidPositions(validPositions);

    unsigned int count = this->fillSystemRec(A,b,waveletMapper,validPositions,0u);

    assert(count == validPositions.size());
}
        
template <typename T>
unsigned int TreeNode<T>::fillSystemRec(Eigen::MatrixXf &A, Eigen::VectorXf &b, 
        const WaveletMapper<T> &waveletMapper, 
        const std::vector<T> &validPositions, 
        unsigned int column) const {
   
    unsigned int newColumn = column;
    unsigned int N = validPositions.size();

    std::cout << "column " << column << std::endl;

    if(this->isValid()) {
        for (unsigned int i = 0; i < N; i++) {
            A(i,column) = (waveletMapper(this->level(), this->offset()))(this->level(), this->offset(), validPositions[i]);
        }

        b(column) = this->value();
        newColumn ++;
    }
    
    TreeNode<T> *child;
    for (unsigned int i = 0; i < _nChilds; i++) {
        child = _childs[i];
        if(child != nullptr) {
            newColumn = child->fillSystemRec(A,b,waveletMapper,validPositions, newColumn);
        }
    }

    return newColumn;
}

#endif /* end of include guard: TREENODE_H */
