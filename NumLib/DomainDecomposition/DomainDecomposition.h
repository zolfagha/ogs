
#pragma once

#include <vector>

#include "Base/BidirectionalMap.h"

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MeshLib/Core/IMesh.h"


namespace NumLib
{

/**
 * \brief Node-based decomposed sub domain
 */
class NodeBasedSubDomain
{
private:
    MathLib::ILinearEquations *_linear_eqs;

    MeshLib::IMesh *_local_msh;
    Base::BidirectionalMap<size_t, size_t> _map_global2localNodeId;
    std::set<size_t> _ghost_nodes;
    std::vector<size_t> _list_dirichlet_bc_id;
    std::vector<double> _list_dirichlet_bc_value;
public:
    void setupEqs()
    {
        size_t n_nodes;
    }
    void assemble() 
    {
        std::vector<int> ele;
        for (size_t j=0; j<ele.size(); j++) {
            // get dof map
            // element assembly
            // add active entries to sub eqs?
        }
        // add to global eqs?
    }
};



class NodeBasedDecomposedDomain
{
private:
    MathLib::ILinearEquations *_linear_eqs;
public:
    size_t getNumberOfSubDomains();
    NodeBasedSubDomain* getSubDomain(size_t i);
    void solve()
    {
        _linear_eqs->solve();
    }
    void doPost()
    {
        // for each sub domain
        for (size_t i_dom=0; i_dom<getNumberOfSubDomains(); i_dom++) {
            // for each elements
            // calculate integration point values
        }

    }

    void assemble() 
    {
        for (size_t i_dom=0; i_dom<getNumberOfSubDomains(); i_dom++) {
            NodeBasedSubDomain* subd = this->getSubDomain(i_dom);
            subd->assemble();
        }
    }


};

class ElementBasedSubDomain
{
private:
    size_t _id;
    MeshLib::IMesh *_msh;
    std::vector<size_t> _nodes_inner;
    std::vector<size_t> _nodes_border;
    std::vector<size_t> _nodes_all;
    std::vector<size_t> _elements;
    std::vector<size_t> _border_nodes;


    void createSparseTable();

public:
    void setMesh(MeshLib::IMesh *msh) {_msh = msh;};
    void createNodes()
    {

    }
    void createElements();
    void createNodeConnectedNodes();
    void createEQS();
    void initializeEQS();
    void configEQS();
    void assemble();
    void applySourceTerns();
    void applyBoundaryConditions();

    double dotBorderVec(double* vec_x, double* vec_y);
    double dotInterior(double*, double*);

    void setLocal2Global(double* local, double* global, size_t n);
    void setGlobal2Local(double* global, double* local, size_t n);
    void setLocal2Border(double* local, double* border);
    void setBorder2Local(double* border, double* local, size_t n);
    void catInnerX(double* global_x, const double* local_x, const long n);
};


class ElementBasedDomainDecompostion
{
private:
    MeshLib::IMesh* _msh;
    std::vector<ElementBasedSubDomain*> _subdomains;

    void findNodesOnInterface();

public:
    void create()
    {
        size_t n_dom = _subdomains.size();
        for (size_t i=0; i<n_dom; i++) {
            ElementBasedSubDomain* sub = _subdomains[i];
            sub->setMesh(_msh);
            sub->createNodes();
            sub->createElements();
        }
        findNodesOnInterface();
        for (size_t i=0; i<n_dom; i++) {
            ElementBasedSubDomain* sub = _subdomains[i];
            sub->createNodeConnectedNodes();
            sub->createEQS();
        }
    }

    void assemble()
    {
        size_t n_dom = _subdomains.size();
        for (size_t i=0; i<n_dom; i++) {
            ElementBasedSubDomain* sub = _subdomains[i];
            sub->initializeEQS();
            sub->assemble();
            sub->applySourceTerns();
            sub->applyBoundaryConditions();
            // add into global
        }
    }
};
}
