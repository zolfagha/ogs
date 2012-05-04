
#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientLinearEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientResidualLocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "Tests/Geo/Material/PorousMedia.h"
#include "Tests/Geo/Material/Solid.h"

namespace Geo
{

class FemLinearElasticLinearLocalAssembler: public NumLib::IElementWiseTransientLinearEQSLocalAssembler
{
private:
	PorousMedia* _pm;
	FemLib::LagrangianFeObjectContainer* _feObjects;
public:
    typedef typename NumLib::IElementWiseTransientLinearEQSLocalAssembler::LocalVectorType LocalVectorType;
    typedef typename NumLib::IElementWiseTransientLinearEQSLocalAssembler::LocalMatrixType LocalMatrixType;

    FemLinearElasticLinearLocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, PorousMedia &pm)
	: _pm(&pm), _feObjects(&feObjects)
	{
	};

	virtual ~FemLinearElasticLinearLocalAssembler() {};

protected:
    virtual void assembly(const NumLib::TimeStep &time,  MeshLib::IElement &e, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalEquationType &eqs)
    {
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

		const size_t dim = e.getDimension();
		const size_t n_strain_components = (dim==2 ? 4 : 6);
		const size_t nnodes = e.getNumberOfNodes();
		Solid *solidphase = _pm->solidphase;

		LocalMatrixType matB(n_strain_components, nnodes*dim);
		LocalMatrixType matD(n_strain_components, n_strain_components);
		LocalMatrixType matDB(matB.getNRows(), matB.getNCols());

		LocalEquationType::MatrixType &localK = *eqs.getA();
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        	q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

            // set D
            matD = .0;
            double nv = solidphase->poisson_ratio;
            double E = solidphase->Youngs_modulus;
            double Lambda, G, K;
            calculateLameConstant(nv, E, Lambda, G, K);
            setElasticConsitutiveTensor(dim, Lambda, G, K, matD);

            // set B
            LocalMatrixType &dN = *fe->getGradBasisFunction();
            matB = .0;
            for (size_t in=0; in<nnodes; in++) {
            	size_t offset = in*dim;
            	if (dim==2) {
                	matB(0,offset) = dN(0,in);
                	matB(1,offset+1) = dN(1,in);
                	matB(3,offset) = dN(0,in);
                	matB(3,offset+1) = dN(1,in);
            	}
            }
        	// K += B^T * D * B
            matD.multiply(matB, matDB);
            matB.transposeAndMultiply(matDB, localK);

            // RHS
        }
    }
};

class FemLinearElasticResidualLocalAssembler : public NumLib::IElementWiseTransientResidualLocalAssembler
{
public:
    virtual ~FemLinearElasticResidualLocalAssembler() {};

    /// assemble a local residual for the given element
    /// @param time			time step
    /// @param e			element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param local_r		local residual
    virtual void assembly(const NumLib::TimeStep &time,  MeshLib::IElement &e, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalVectorType &local_r)
    {

    }
};

class FemLinearElasticJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
private:
	PorousMedia* _pm;
	FemLib::LagrangianFeObjectContainer* _feObjects;
public:
	FemLinearElasticJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, PorousMedia &pm)
	: _pm(&pm), _feObjects(&feObjects)
	{
	};

	void assembly(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/,  LocalMatrixType &localJ)
	{
		FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        	q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

        	double k;
        	_pm->hydraulic_conductivity->eval(real_x, k);

    		//fe->integrateWxN(_pm->storage, localM);
    		fe->integrateDWxDN(j, k, localJ);
        }
	}
};


} //end
