/* File : example.c */

#include "example.h"

#include "Tests/TestExamples.h"
#include "MeshLib/Tools/MeshGenerator.h"

#define M_PI 3.14159265358979323846

/* Move the shape to a new location */
void Shape::move(double dx, double dy) {
  x += dx;
  y += dy;
}

int Shape::nshapes = 0;

double Circle::area(void) {
  return M_PI*radius*radius;
}

double Circle::perimeter(void) {
  return 2*M_PI*radius;
}

double Square::area(void) {
  return width*width;
}

double Square::perimeter(void) {
  return 4*width;
}

void gw1(double len, size_t div, std::vector<double> *x, std::vector<double> *results)
{
    GWFemTest gw;
    MeshLib::IMesh *msh = MeshLib::MeshGenerator::generateRegularQuadMesh(len, div, .0, .0, .0);
    gw.define(msh);
    //#Solve
    GWFemTest::calculateHead(gw);

    DiscreteLib::DiscreteVector<double>* h = gw.head->getNodalValues();

    x->resize(msh->getNumberOfNodes());
    for (size_t i=0; i<x->size(); i++)
        (*x)[i] = msh->getNodeCoordinatesRef(i)->getData()[0];

    results->resize(h->size());
    for (size_t i=0; i<results->size(); i++)
        (*results)[i] = (*h)[i];
}

void range(std::vector<int> *rangevec)
{
    int i;
    for (i=0; i< rangevec->size(); i++)
        (*rangevec)[i] = i;
}
