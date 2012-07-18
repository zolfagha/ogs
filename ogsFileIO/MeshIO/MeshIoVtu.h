
#pragma once

#include <string>
#include <vector>

#include "MeshLib/Core/IMesh.h"


class MeshIoVtu
{
public:
    enum VTK_XML_DATA_TYPE { Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32,
                         Float64 };

    explicit MeshIoVtu(bool binary_mode);
    virtual ~MeshIoVtu(void){}

public:
    bool WriteXMLUnstructuredGrid(const std::string &vtkfile,
                                  MeshLib::IMesh &msh);

protected:
    void InitializeVTU();

    unsigned char GetVTKCellType(const MeshLib::ElementShape::type ele_type);

    bool WriteDataArrayHeader(std::fstream &fin,
                              VTK_XML_DATA_TYPE data_type,
                              const std::string &str_name,
                              int nr_components,
                              const std::string &str_format,
                              long offset = -1);

    bool WriteDataArrayFooter(std::fstream &fin);

    bool WriteMeshNodes(std::fstream &fin,
                               bool output_data,
                               MeshLib::IMesh& m_msh,
                               long &offset);

    bool WriteMeshElementConnectivity(std::fstream &fin,
                                             bool output_data,
                                             MeshLib::IMesh& m_msh,
                                             long &offset,
                                             long &sum_ele_components);

    bool WriteMeshElementOffset(std::fstream &fin,
                                       bool output_data,
                                       MeshLib::IMesh& m_msh,
                                       long &offset);

    bool WriteMeshElementType(std::fstream &fin,
                                     bool output_data,
                                     MeshLib::IMesh& m_msh,
                                     long &offset);

    bool WriteElementGroupID(std::fstream &fin,
                                  bool output_data,
                                  MeshLib::IMesh& m_msh,
                                  long &offset);

protected:
    //for binary output
    bool useBinary;
    bool isInitialized;
    VTK_XML_DATA_TYPE type_UChar;
    VTK_XML_DATA_TYPE type_Int;
    VTK_XML_DATA_TYPE type_UInt;
    VTK_XML_DATA_TYPE type_Long;
    VTK_XML_DATA_TYPE type_Double;
    int SIZE_OF_BLOCK_LENGTH_TAG;
    bool isLittleEndian;                  //Endian(byte order)

};
