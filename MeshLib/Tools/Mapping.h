
#pragma once

#include <vector>

#include "Base/MemoryTools.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/MathTools.h"
#include "GeoLib/Core/Point.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElement.h"
#include "MeshLib/Core/CoordinateSystem.h"

namespace MeshLib
{

/**
 * \brief Interface of methods mapping coordinates of elements.
 *
 *
 */
class IElementCoordinatesMapping
{
public:
    /// get mapped coordinates of nodes
    virtual GeoLib::Point* getNodePoint(size_t node_id) = 0;
};

/**
 * \brief EleMapInvariant keep original coordinates.
 */
class EleMapInvariant : public IElementCoordinatesMapping
{
public:
    EleMapInvariant(const IMesh* msh, IElement* e) 
    {
        _msh = msh;
        _e = e;
    };

    GeoLib::Point* getNodePoint(size_t local_id) 
    {
        return (GeoLib::Point*)_msh->getNodeCoordinatesRef(_e->getNodeID(local_id));
    }

private:
    const IMesh* _msh;
    IElement *_e;
};

/**
 * \brief Mapping local coordinates of elements.
 */
class EleMapLocalCoordinates : public IElementCoordinatesMapping 
{
public:
    ///
    EleMapLocalCoordinates(const IMesh* msh, IElement* e, const CoordinateSystem* coordinate_system) : _matR(0)
    {
        assert (e->getDimension() <= coordinate_system->getDimension());
        _msh = msh;

        if (e->getDimension()==coordinate_system->getDimension()) {
            flip(e, coordinate_system);
        } else if (e->getDimension() < coordinate_system->getDimension()) {
            rotate(e, coordinate_system);
        }
    };

    virtual ~EleMapLocalCoordinates()
    {
        Base::destroyStdVectorWithPointers(_point_vec);
        if (_matR!=0)
            delete _matR;
    }

    virtual GeoLib::Point* getNodePoint(size_t node_id) {
        return _point_vec[node_id];
    }

private:
    const IMesh* _msh;
    std::vector<GeoLib::Point*> _point_vec;
    MathLib::Matrix<double> *_matR;

    ///
    void flip(IElement* e, const CoordinateSystem* coordinate_system)
    {
        switch(coordinate_system->getType())
        {
        case CoordinateSystem::Y:
            {
                assert(e->getDimension()==1);
                for(size_t i = 0; i < e->getNumberOfNodes(); i++)
                {
                    const GeoLib::Point *p = _msh->getNodeCoordinatesRef(e->getNodeID(i));
                    _point_vec.push_back(new GeoLib::Point((*p)[1], (*p)[0], (*p)[2]));
                }
            }
            break;
        case CoordinateSystem::Z:
            {
                assert(e->getDimension()==1);
                for(size_t i = 0; i < e->getNumberOfNodes(); i++)
                {
                    const GeoLib::Point *p = _msh->getNodeCoordinatesRef(e->getNodeID(i));
                    _point_vec.push_back(new GeoLib::Point((*p)[2], (*p)[1], (*p)[0]));
                }
            }
            break;
        case CoordinateSystem::XZ:
            {
                assert(e->getDimension()==2);
                for(size_t i = 0; i < e->getNumberOfNodes(); i++)
                {
                    const GeoLib::Point *p = _msh->getNodeCoordinatesRef(e->getNodeID(i));
                    _point_vec.push_back(new GeoLib::Point((*p)[0], (*p)[2], (*p)[1]));
                }
            }
            break;
        default:
            {
                for(size_t i = 0; i < e->getNumberOfNodes(); i++)
                {
                    const GeoLib::Point *p = _msh->getNodeCoordinatesRef(e->getNodeID(i));
                    _point_vec.push_back(new GeoLib::Point(p->getData()));
                }
            }
        }
    }

    ///
    void rotate(IElement* e, const CoordinateSystem* coordinate_system)
    {
        _point_vec.resize(e->getNumberOfNodes());

        std::vector<size_t> vec_node_id;
        e->getNodeIDList(vec_node_id);
        std::vector<GeoLib::Point> vec_pt;
        _msh->getListOfNodeCoordinates(vec_node_id, vec_pt);

        if (_matR==0)
            getRotationMatrix(e, coordinate_system, vec_pt);

        double const* const coords_node_0 (vec_pt[0].getData());
        double dx[3];
        for(size_t i = 0; i < e->getNumberOfNodes(); i++)
        {
            double const* const coords_node_i (vec_pt[i].getData());
            dx[0] = (coords_node_i[0] - coords_node_0[0]);
            dx[1] = (coords_node_i[1] - coords_node_0[1]);
            dx[2] = (coords_node_i[2] - coords_node_0[2]);

            GeoLib::Point *p = new GeoLib::Point();
            _matR->axpy(1.0, dx, 0.0, (double*)p->getData());
            _point_vec[i] = p;
        }
    };

    void getRotationMatrix(const IElement* e, const CoordinateSystem* coordinate_system, const std::vector<GeoLib::Point> &vec_pt)
    {
        double xx[3];
        double yy[3];
        double zz[3];
        _matR = new MathLib::Matrix<double>(3,3);
        if (e->getDimension() == 1)
        {
            // x"_vec
            double const* const pnt0(vec_pt[0].getData());
            double const* const pnt1(vec_pt[1].getData());
            xx[0] = pnt1[0] - pnt0[0];
            xx[1] = pnt1[1] - pnt0[1];
            xx[2] = pnt1[2] - pnt0[2];
            MathLib::normalizeVector(xx, 3);
            // an arbitrary vector
            for (size_t i = 0; i < 3; i++)
                yy[i] = 0.0;
            if (fabs(xx[0]) > 0.0 && fabs(xx[1]) + fabs(xx[2]) < DBL_MIN)
                yy[2] = 1.0;
            else if (fabs(xx[1]) > 0.0 && fabs(xx[0]) + fabs(xx[2]) < DBL_MIN)
                yy[0] = 1.0;
            else if (fabs(xx[2]) > 0.0 && fabs(xx[0]) + fabs(xx[1]) < DBL_MIN)
                yy[1] = 1.0;
            else
            {
                for (size_t i = 0; i < 3; i++)
                    if (fabs(xx[i]) > 0.0)
                    {
                        yy[i] = -xx[i];
                        break;
                    }
            }
            // z"_vec
            MathLib::crossProd(xx, yy, zz);
            MathLib::normalizeVector(zz, 3);
            // y"_vec
            MathLib::crossProd(zz, xx, yy);
            MathLib::normalizeVector(yy, 3);
        }
        else if (e->getDimension()==2)
        {
            // x"_vec
            //			xx[0] = nodes[1]->X() - nodes[0]->X();
            //			xx[1] = nodes[1]->Y() - nodes[0]->Y();
            //			xx[2] = nodes[1]->Z() - nodes[0]->Z();
            double const* const pnt0(vec_pt[0].getData());
            double const* const pnt1(vec_pt[1].getData());
            xx[0] = pnt1[0] - pnt0[0];
            xx[1] = pnt1[1] - pnt0[1];
            xx[2] = pnt1[2] - pnt0[2];
            MathLib::normalizeVector(xx, 3);
            // a vector on the plane
            //			yy[0] = nodes[2]->X() - nodes[1]->X();
            //			yy[1] = nodes[2]->Y() - nodes[1]->Y();
            //			yy[2] = nodes[2]->Z() - nodes[1]->Z();
            double const* const pnt2(vec_pt[2].getData());
            yy[0] = pnt2[0] - pnt1[0];
            yy[1] = pnt2[1] - pnt1[1];
            yy[2] = pnt2[2] - pnt1[2];
            // z"_vec. off plane
            MathLib::crossProd(xx, yy, zz);
            MathLib::normalizeVector(zz, 3);
            // y"_vec
            MathLib::crossProd(zz, xx, yy);
            MathLib::normalizeVector(yy, 3);
        }
        for (size_t i = 0; i < 3; i++)
        {
            (*_matR)(i, 0) = xx[i];
            (*_matR)(i, 1) = yy[i];
            (*_matR)(i, 2) = zz[i];
        }
    }
};

}
