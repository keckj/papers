
#ifndef WAVELETTREE_H
#define WAVELETTREE_H

#include <Eigen/Dense>
#include "treeNode.hpp"
#include "waveletMapper.hpp"

template <typename T>
class WaveletTree {

    public:
        explicit WaveletTree(const TreeNode<T> *tree, const Wavelet<T> &wavelet, T coef);
        virtual ~WaveletTree();

        int level() const;
        int offset() const;
        
        bool isValid() const;

        unsigned int nChilds() const; 
        WaveletTree<T>* getChild(unsigned int childId) const;
        WaveletTree<T>*& child(unsigned int childId);

        T operator()(T x) const;
        
        void plot(Gnuplot &gp) const;
        
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
    _isValid(tree->isValid()), _coefficient(coef), 
    _wavelet(wavelet) {
}
       
template <typename T>
WaveletTree<T>* WaveletTree<T>::makeTree(const TreeNode<T> *node, 
        const WaveletMapper<T> &mapper, 
        const Eigen::VectorXf &coefficents) {
    
    unsigned int count = 0u;
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
    
    const TreeNode<T> *child;
    for (unsigned int i = 0; i < node->nChilds(); i++) {
        child = node->getChild(i);
        if(child != nullptr)
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
void WaveletTree<T>::plot(Gnuplot &gp) const {
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
void WaveletTree<T>::plotRec(std::vector<std::tuple<T,T>> &pts, unsigned int jmax) const {
   
    std::tuple<T,T> c = std::make_tuple<T,T>(this->position(), jmax + 1u - this->level());

    //insert child points
    WaveletTree<T> *node; 
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
T WaveletTree<T>::operator()(T x) const {
}

#endif /* end of include guard: WAVELETTREE_H */
