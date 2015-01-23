
#ifndef WAVELETTREE_H
#define WAVELETTREE_H

#include <Eigen/Dense>
#include "treeNode.hpp"
#include "plotBox.hpp"
#include "waveletMapper.hpp"

template <typename T>
class WaveletTree {

    public:
        explicit WaveletTree(const TreeNode<T> *tree, const Wavelet<T> &wavelet, T coef);
        virtual ~WaveletTree();

        const Wavelet<T>& wavelet() const;

        int level() const;
        int offset() const;
        
        bool isValid() const;

        Interval<T> support() const;
        T inf() const;
        T sup() const;

        T coefficient() const;

        unsigned int nChilds() const; 
        WaveletTree<T>* getChild(unsigned int childId) const;
        WaveletTree<T>*& child(unsigned int childId);

        T operator()(T x) const;
        
        void plot(Gnuplot &gp, const PlotBox<T> &box, unsigned int nPoints) const;
        
        static  WaveletTree<T>* makeTree(const TreeNode<T> *node, 
                const WaveletMapper<T> &mapper, 
                const Eigen::VectorXf &coefficents);

        
    protected:
        const unsigned int _level;
        const int _offset;

        const unsigned int _nChilds;
        WaveletTree<T> **_childs;
    
        const bool _isValid;
        const T _coefficient;
        const Wavelet<T> &_wavelet;
        
        void plotRec(std::vector<std::tuple<T,T>> &pts, unsigned int jmax) const;
        
        static WaveletTree<T>* makeTreeRec(const TreeNode<T> *node, 
                const WaveletMapper<T> &mapper, 
                const Eigen::VectorXf &coefficents,
                unsigned int &count);
};


template <typename T>
WaveletTree<T>::WaveletTree(const TreeNode<T> *tree, const Wavelet<T> &wavelet, T coef) :
    _level(tree->level()), _offset(tree->offset()),
    _nChilds(tree->nChilds()), _childs(nullptr),
    _isValid(tree->isValid()&&tree->isFilled()), _coefficient(coef), 
    _wavelet(wavelet) {

    _childs = new WaveletTree*[_nChilds]; 
    for (unsigned int i = 0; i < _nChilds; i++) {
        _childs[i] = nullptr;
    }
}
       
template <typename T>
WaveletTree<T>* WaveletTree<T>::makeTree(const TreeNode<T> *node, 
        const WaveletMapper<T> &mapper, 
        const Eigen::VectorXf &coefficents) {
    
        unsigned int count = 0u;
        std::cout << node->interval() << std::endl;
        return WaveletTree<T>::makeTreeRec(node, mapper, coefficents, count);
}

template <typename T>
WaveletTree<T>* WaveletTree<T>::makeTreeRec(const TreeNode<T> *node, 
        const WaveletMapper<T> &mapper, 
        const Eigen::VectorXf &coefficents,
        unsigned int &count) {
   
    WaveletTree<T>* waveNode;

    if(node->isValid()) {
        waveNode = new WaveletTree<T>(node, mapper(node->level(), node->offset()),coefficents[count]);
        count++;
    }
    else {
        waveNode = new WaveletTree<T>(node, mapper(node->level(), node->offset()),T(0));
    }

    if(node->isWeakLinked())
        return waveNode;
    
    const TreeNode<T> *child;
    for (unsigned int i = 0; i < node->nChilds(); i++) {
        child = node->getChild(i);
        if(child != nullptr && child->isFilled())
            waveNode->child(i) = makeTreeRec(child, mapper, coefficents, count);
    }

    return waveNode;
}

template <typename T>
WaveletTree<T>::~WaveletTree() {
    for (unsigned int i = 0; i < _nChilds; i++) {
        delete _childs[i];
    }
    delete [] _childs;
}

template <typename T>
 int WaveletTree<T>::level() const {
	 return _level;
}

template <typename T>
 int WaveletTree<T>::offset() const {
	 return _offset;
}

template <typename T>
 bool WaveletTree<T>::isValid() const {
	 return _isValid;
}

template <typename T>
unsigned int WaveletTree<T>::nChilds() const {
	 return _nChilds;
}

template <typename T>
WaveletTree<T>* WaveletTree<T>::getChild(unsigned int childId) const {
    return _childs[childId];
}
        
template <typename T>
WaveletTree<T>*& WaveletTree<T>::child(unsigned int childId) {
    return _childs[childId];
}

template <typename T>
void WaveletTree<T>::plot(Gnuplot &gp, const PlotBox<T> &box, unsigned int nPoints) const {

    std::vector<std::tuple<T,T>> pts;

    Interval<T> interval(box.xmin(), box.xmax());

    T dx = interval.length()/(nPoints - 1u);
    for (unsigned int i = 0; i < nPoints; i++) {
        T x = box.xmin() + i*dx;
        T y = this->operator()(x);
        pts.push_back(std::make_tuple(x, y));
    }
    
    gp << box;
    gp << "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5\n";
    gp << "plot '-' with lines ls 1 title 'wavelet'\n";
    gp.send1d(pts);
}

template <typename T>
T WaveletTree<T>::operator()(T x) const {

   //std::cout << "level " << this->level() << "\toffset " << this->offset() << " support " << this->support() << std::endl;

   if(!this->support().contains(x) && this->level() >= 0)
           return T(0);

   T res;
   if(this->isValid()) {
       res = this->coefficient()*this->wavelet()(this->level(), this->offset(), x);
    }
   else { 
       res = T(0); 
    }
    
   //compute child values
    WaveletTree<T> *child; 
    for (unsigned int i = 0; i < this->nChilds(); i++) {
        child = _childs[i];
        if(child != nullptr) {
            res += child->operator()(x);
        }
    }

    return res;
}
        
template <typename T>
Interval<T> WaveletTree<T>::support() const {
    return this->wavelet().support(this->level(), this->offset());
}

template <typename T>
T WaveletTree<T>::inf() const {
    return this->support().inf;
}

template <typename T>
T WaveletTree<T>::sup() const {
    return this->support().sup;
}
        
template <typename T>
const Wavelet<T>& WaveletTree<T>::wavelet() const {
    return _wavelet;
}
        
template <typename T>
T WaveletTree<T>::coefficient() const {
    return _coefficient;
}

#endif /* end of include guard: WAVELETTREE_H */
