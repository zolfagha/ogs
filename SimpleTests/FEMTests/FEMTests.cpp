
#include <iostream>
//#include <cmath>
//#include "Point.h"
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
    std::string strMeshFile = "";
    std::cout << "Input a mesh file name:" << endl;
    cin >> strMeshFile;
    vector<IMesh*> vec_mesh; 
    MeshIOOGS::readMesh(strMeshFile, vec_mesh);
    //MeshGenerator::generateRegularMesh(2, 100, 100, 0.0, 0.0, 0.0, msh);
    //MeshIOLegacyVtk::WriteAsciiFile("test.vtk", msh);
    if (vec_mesh.size()==0) {
      std::cout << "Fail to read a mesh file: " << strMeshFile << std::endl;
      return 0;
    }

    std::cout << "== Mesh ==" << endl;
    for (size_t i=0; i<vec_mesh.size(); i++) {
      std::cout << "#Mesh " << i+1 << endl;
      std::cout << "Number of nodes   : " << vec_mesh[i]->getNumberOfNodes() << endl;
      std::cout << "Number of elements: " << vec_mesh[i]->getNumberOfElements() << endl;
    }
    std::cout << endl;



    //prepare EQS, DoF

    //assembly EQS
    for (size_t i=0; i<0; i++) {

    }

    //apply BC/ST

    //solve EQS

    //

    return 0;
}

#endif
