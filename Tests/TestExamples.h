
#pragma once

#include <vector>

#include "Base/BidirectionalMap.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Assembler/ElementLocalAssembler.h"

struct DiscreteExample1
{
    std::vector<size_t> list_dirichlet_bc_id;
    std::vector<double> list_dirichlet_bc_value;
    static const size_t dim_eqs = 9;
    std::vector<double> exH;

    DiscreteExample1()
    {
        size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        list_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        list_dirichlet_bc_value.resize(6);
        fill(list_dirichlet_bc_value.begin(), list_dirichlet_bc_value.begin()+3, .0);
        fill(list_dirichlet_bc_value.begin()+3, list_dirichlet_bc_value.end(), 1.0);
        exH.resize(9);
        for (size_t i=0; i<9; i++) {
            if (i%3==0) exH[i] = 1.0;
            if (i%3==1) exH[i] = 0.5;
            if (i%3==2) exH[i] = 0.;
        }
    }

    void setLocalDirichletBC(const Base::BidirectionalMap<size_t, size_t> &map_global2localNodeId, std::vector<size_t> &local_dirichlet_bc_id, std::vector<double> &local_dirichlet_bc_value)
    {
        for (size_t i=0; i<list_dirichlet_bc_id.size(); i++) {
            if (map_global2localNodeId.countInA(list_dirichlet_bc_id[i])>0) {
                size_t local_id = map_global2localNodeId.mapAtoB(list_dirichlet_bc_id[i]);
                local_dirichlet_bc_id.push_back(local_id);
                local_dirichlet_bc_value.push_back(list_dirichlet_bc_value[i]);
            }
        }
    }


    class TestElementAssembler : public DiscreteLib::IElemenetLocalAssembler
    {
        MathLib::Matrix<double> _m;
    public:
        TestElementAssembler()
        {
            _m.resize(4,4);
            _m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0; 
            _m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
            _m(2,2) = 4.0; _m(2,3) = -1.0;
            _m(3,3) = 4.0;
            for (size_t i=0; i<4; i++)
                for (size_t j=0; j<i; j++) _m(i,j) = _m(j,i);
            _m *= 1.e-11/6.0;
        }
        void assembly(MeshLib::IElement &e, MathLib::DenseLinearEquations &eqs)
        {
            (*eqs.getA()) = _m;
        }
    };
};
