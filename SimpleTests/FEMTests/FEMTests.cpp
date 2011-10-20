
#include <iostream>
#include <cmath>
#include "Point.h"
#include "Mesh.h"
#include "MeshIO.h"
#include "MeshGenerator.h"

using namespace std;
using namespace MeshLib;

#define TEST2

#ifdef TEST1
int main(int argc, char *argv[])
{
    cout << "Hello World!" << endl;
    StaticUnstructuredMesh<3,10,10> mesh;
    cout << "Nr. Nodes=" << mesh.getNumberOfNodes() << endl;
    cout << "Nr. Elements=" << mesh.getNumberOfElements() << endl;
    
	return 0;
}

#elif defined TEST2
int main(int argc, char *argv[])
{
    //set IC,BC,ST
    //create mesh
    UnstructuredMesh msh;
    MeshGenerator::generateRegularMesh(2, 100, 100, 0.0, 0.0, 0.0, msh);

    MeshIOLegacyVtk::WriteAsciiFile("test.vtk", msh);

    //prepare EQS

    //assembly EQS
    for (size_t i=0; i<0; i++) {

    }

    //apply BC/ST

    //solve EQS

    //

    return 0;
}

#endif
