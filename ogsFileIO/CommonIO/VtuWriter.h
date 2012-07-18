
#pragma once

#include <string>
#include <vector>

#include "MeshLib/Core/IMesh.h"
#include "NumLib/Function/TXFunction.h"

class VtuWriter
{
public:
    enum VTK_XML_DATA_TYPE { Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32,
                         Float64 };

    typedef std::pair<std::string, NumLib::ITXFunction*> PointData;
    typedef std::pair<std::string, NumLib::ITXFunction*> CellData;

    explicit VtuWriter(bool binary_mode);
    virtual ~VtuWriter(void){}

public:
    bool write(const std::string &vtkfile,
                                  MeshLib::IMesh &msh, std::vector<PointData> &nodal, std::vector<CellData> &elemental);

protected:
    void initialize();

    unsigned char getVTKCellType(const MeshLib::ElementShape::type ele_type);

    bool writeDataArrayHeader(std::fstream &fin,
                              VTK_XML_DATA_TYPE data_type,
                              const std::string &str_name,
                              int nr_components,
                              const std::string &str_format,
                              long offset = -1);

    bool writeDataArrayFooter(std::fstream &fin);

    bool writePoints(std::fstream &fin,
                               bool output_data,
                               MeshLib::IMesh& m_msh,
                               long &offset);

    bool writeCellConnectivity(std::fstream &fin,
                                             bool output_data,
                                             MeshLib::IMesh& m_msh,
                                             long &offset,
                                             long &sum_ele_components);

    bool writeCellOffset(std::fstream &fin,
                                       bool output_data,
                                       MeshLib::IMesh& m_msh,
                                       long &offset);

    bool writeCellType(std::fstream &fin,
                                     bool output_data,
                                     MeshLib::IMesh& m_msh,
                                     long &offset);
    
    bool writePointData(std::fstream &fin,
                                bool output_data,
                                std::vector<PointData> &nodal,
                                MeshLib::IMesh &m_msh,
                                long &offset);
    bool writeCellData(std::fstream &fin,
                                  bool output_data,
                                  std::vector<CellData> &elemental,
                                    MeshLib::IMesh &m_msh,
                                  long &offset);

protected:
    //for binary output
    bool _useBinary;
    bool _isInitialized;
    //Endian(byte order)
    bool _isLittleEndian;
    
    VTK_XML_DATA_TYPE type_UChar;
    VTK_XML_DATA_TYPE type_Int;
    VTK_XML_DATA_TYPE type_UInt;
    VTK_XML_DATA_TYPE type_Long;
    VTK_XML_DATA_TYPE type_Double;
    int SIZE_OF_BLOCK_LENGTH_TAG;
    
};
