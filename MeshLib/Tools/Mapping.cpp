
#include "Mapping.h"

#include <limits>
#include "MathLib/MathTools.h"

namespace MeshLib
{

///
void EleMapLocalCoordinates::flip(IElement &ele, const CoordinateSystem &coordinate_system)
{
    IElement* e = &ele;
    switch(coordinate_system.getType())
    {
    case CoordinateSystemType::Y:
        {
            assert(e->getDimension()==1);
            for(size_t i = 0; i < e->getNumberOfNodes(); i++)
            {
                const GeoLib::Point *p = _msh->getNodeCoordinatesRef(e->getNodeID(i));
                _point_vec.push_back(new GeoLib::Point((*p)[1], (*p)[0], (*p)[2]));
            }
        }
        break;
    case CoordinateSystemType::Z:
        {
            assert(e->getDimension()==1);
            for(size_t i = 0; i < e->getNumberOfNodes(); i++)
            {
                const GeoLib::Point *p = _msh->getNodeCoordinatesRef(e->getNodeID(i));
                _point_vec.push_back(new GeoLib::Point((*p)[2], (*p)[1], (*p)[0]));
            }
        }
        break;
    case CoordinateSystemType::XZ:
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
void EleMapLocalCoordinates::rotate(IElement &ele, const CoordinateSystem &coordinate_system)
{
    IElement* e = &ele;
    _point_vec.resize(e->getNumberOfNodes());

    std::vector<size_t> vec_node_id;
    e->getNodeIDList(vec_node_id);
    std::vector<GeoLib::Point> vec_pt;
    _msh->getListOfNodeCoordinates(vec_node_id, vec_pt);

    if (_matR2original==0) {
        getRotationMatrixToOriginal(*e, coordinate_system, vec_pt);
        _matR2local = new MathLib::Matrix<double>(3,3);
        _matR2original->transpose(*_matR2local);
    }

    double const* const coords_node_0 (vec_pt[0].getData());
    double dx[3];
    for(size_t i = 0; i < e->getNumberOfNodes(); i++)
    {
        double const* const coords_node_i (vec_pt[i].getData());
        dx[0] = (coords_node_i[0] - coords_node_0[0]);
        dx[1] = (coords_node_i[1] - coords_node_0[1]);
        dx[2] = (coords_node_i[2] - coords_node_0[2]);

        GeoLib::Point *p = new GeoLib::Point();
        _matR2local->axpy(1.0, dx, 0.0, (double*)p->getData());
        _point_vec[i] = p;
    }
};

// x=Rx' where x is original coordinates and x' is local coordinates
void EleMapLocalCoordinates::getRotationMatrixToOriginal(const IElement &ele, const CoordinateSystem &coordinate_system, const std::vector<GeoLib::Point> &vec_pt)
{
    const IElement* e = &ele;
    double xx[3];
    double yy[3];
    double zz[3];
    assert(_matR2original==0);
    _matR2original = new MathLib::Matrix<double>(3,3);

    if (e->getDimension() == 1) {
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

        if (fabs(xx[0]) > 0.0 && fabs(xx[1]) + fabs(xx[2]) < std::numeric_limits<double>::epsilon()) {
            yy[2] = 1.0;
        } else if (fabs(xx[1]) > 0.0 && fabs(xx[0]) + fabs(xx[2]) < std::numeric_limits<double>::epsilon()) {
            yy[0] = 1.0;
        } else if (fabs(xx[2]) > 0.0 && fabs(xx[0]) + fabs(xx[1]) < std::numeric_limits<double>::epsilon()) {
            yy[1] = 1.0;
        } else {
            for (size_t i = 0; i < 3; i++) {
                if (fabs(xx[i]) > 0.0) {
                    yy[i] = -xx[i];
                    break;
                }
            }
        }
        // z"_vec
        MathLib::crossProd(xx, yy, zz);
        MathLib::normalizeVector(zz, 3);
        // y"_vec
        MathLib::crossProd(zz, xx, yy);
        MathLib::normalizeVector(yy, 3);

    } else if (e->getDimension()==2) {
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

    for (size_t i=0; i<3; ++i) {
        (*_matR2original)(i, 0) = xx[i];
        (*_matR2original)(i, 1) = yy[i];
        (*_matR2original)(i, 2) = zz[i];
    }
}

}
