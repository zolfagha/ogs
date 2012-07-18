
#include "VtuWriter.h"

#include <fstream>

#include "BaseLib/SystemTools.h"
#include "BaseLib/FileTools.h"
#include "MeshLib/Core/CoordinateSystem.h"
#include "MeshLib/Core/MeshGeometricProperties.h"

using namespace MeshLib;

VtuWriter::VtuWriter(bool binary_mode)
: _useBinary(binary_mode), _isInitialized(false), _isLittleEndian(true), SIZE_OF_BLOCK_LENGTH_TAG(0)
{
    this->initialize();
}

bool VtuWriter::write(const std::string &vtkfile,
                                    MeshLib::IMesh &msh, std::vector<PointData> &nodal_values, std::vector<CellData> &elemental_values)
{
    //-------------------------------------------------------------------------
    //# Setup file stream
    //-------------------------------------------------------------------------
    std::fstream fin;
    if (_useBinary)
        fin.open(vtkfile.data(), std::ios::out | std::ios::binary);
    else
        fin.open(vtkfile.data(), std::ios::out);

    if (!fin.good())
    {
        std::cout << "***Warning: Cannot open the output file, " << vtkfile << std::endl;
        return false;
    }

    if (!_useBinary)
    {
        fin.setf(std::ios::scientific,std::ios::floatfield);
        fin.precision(12);
    }

    //-------------------------------------------------------------------------
    //# Output
    //-------------------------------------------------------------------------
    long offset = 0;

    std::string str_format;
    if (!_useBinary)
        str_format = "ascii";
    else
        str_format = "appended";
    bool data_out = !_useBinary;

    //# Header
    fin << "<?xml version=\"1.0\"?>" << std::endl;
    fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
    if (!_useBinary || _isLittleEndian)
        fin << " byte_order=\"LittleEndian\"";
    else
        fin << " byte_order=\"BigEndian\"";
    fin << ">" << std::endl;
    //  fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"  << endl;

    //# Unstructured Grid information
    fin << "  <UnstructuredGrid>" << std::endl;
    fin << "    <Piece NumberOfPoints=\"" << msh.getNumberOfNodes() <<
    "\" NumberOfCells=\"" << msh.getNumberOfElements() << "\">" << std::endl;
    //....................................................................
    // Nodes
    fin << "      <Points>" << std::endl;
    writeDataArrayHeader(fin, type_Double, "", 3, str_format, offset);
    writePoints(fin, data_out, msh, offset);
    writeDataArrayFooter(fin);
    fin << "      </Points>" << std::endl;
    //....................................................................
    // Elements
    fin << "      <Cells>" << std::endl;
    //connectivity
    writeDataArrayHeader(fin, type_Long, "connectivity", 0, str_format, offset);
    long sum_ele_components = 0;
    writeCellConnectivity(fin, data_out, msh, offset, sum_ele_components);
    writeDataArrayFooter(fin);
    //offset
    writeDataArrayHeader(fin, type_Long, "offsets", 0, str_format, offset);
    writeCellOffset(fin, data_out, msh, offset);
    writeDataArrayFooter(fin);
    //type
    writeDataArrayHeader(fin, type_UChar, "types", 0, str_format, offset);
    writeCellType(fin, data_out, msh, offset);
    writeDataArrayFooter(fin);
    fin << "      </Cells>" << std::endl;
    //....................................................................
    // Nodal values
    if (nodal_values.size() > 0)
        fin << "      <PointData Scalars=\"" << nodal_values[0].first << "\">" << std::endl;
    else
        fin << "      <PointData Scalars=\"scalars\">" << std::endl;
    writePointData(fin, data_out, nodal_values, msh, offset);
    fin << "      </PointData>" << std::endl;
    //....................................................................
    // Element values
    fin << "      <CellData>" << std::endl;
    writeCellData(fin, data_out, elemental_values, msh, offset);
    fin << "      </CellData>" << std::endl;
    fin << "    </Piece>" << std::endl;
    fin << "  </UnstructuredGrid>" << std::endl;

    //======================================================================
    // Raw data (for binary mode)
    if (_useBinary)
    {
        fin << "  <AppendedData encoding=\"raw\">" << std::endl;
        fin << "    _";

        //Node
        this->writePoints(fin, true, msh, offset);
        //Element
        //connectivity
        this->writeCellConnectivity(fin, true, msh, offset, sum_ele_components);
        //offset
        this->writeCellOffset(fin, true, msh, offset);
        //type
        this->writeCellType(fin, true, msh, offset);
        // Nodal values
        this->writePointData(fin, true, nodal_values, msh, offset);
        // Elemental values
        this->writeCellData(fin, true, elemental_values, msh, offset);

        fin << std::endl;
        fin << "  </AppendedData>" << std::endl;
    }

    fin << "</VTKFile>" << std::endl;
    fin.close();

    return true;
}

bool VtuWriter::writePointData(std::fstream &fin,
                           bool output_data,
                           std::vector<PointData> &nod_values,
                           IMesh &msh,
                           long &offset)
{
    std::string str_format;
    if (!this->_useBinary)
        str_format = "ascii";
    else
        str_format = "appended";

    //Nodal values
    for (size_t i=0; i<nod_values.size(); i++)
    {
        const std::string &var_name = nod_values[i].first;

        if (!_useBinary || !output_data)
            writeDataArrayHeader(fin, type_Double, var_name, 0, str_format, offset);

        if (output_data)
        {
            if (!_useBinary) {
                fin << "          ";
            } else {
                BaseLib::write_value_binary<unsigned int> (fin, sizeof(double)* msh.getNumberOfNodes());
            }
            for (size_t j = 0; j < msh.getNumberOfNodes(); j++) {
                
                double v = .0;
                nod_values[i].second->eval(NumLib::TXPosition(j, msh.getNodeCoordinatesRef(j)), v);
                if (!_useBinary) {
                    fin << v << " ";
                } else {
                    BaseLib::write_value_binary(fin, v);
                }
            }
            if (!_useBinary) {
                fin << std::endl;
            }
        }
        else
        {
            offset += msh.getNumberOfNodes() * sizeof(double)
                      + SIZE_OF_BLOCK_LENGTH_TAG;
        }

        if (!_useBinary || !output_data)
            writeDataArrayFooter(fin);
    }


    return true;
}

bool VtuWriter::writeCellData(std::fstream &fin,
                             bool output_data,
                             std::vector<CellData> &ele_values,
                             IMesh &msh,
                             long &offset)
{
    std::string str_format;
    if (!this->_useBinary)
        str_format = "ascii";
    else
        str_format = "appended";

    //Element values
    for (int i = 0; i < (int) ele_values.size(); i++)
    {
        if (!_useBinary || !output_data)
            writeDataArrayHeader(fin, this->type_Double,
                    ele_values[i].first, 0, str_format, offset);
        if (output_data)
        {
            if (!_useBinary)
            {
                fin << "          ";
                for (long j = 0; j < (long) msh.getNumberOfElements(); j++) {
                    double v = .0;
                    ele_values[i].second->eval(NumLib::TXPosition(j), v);
                    fin << v << " ";
                }
                fin << std::endl;
            }
            else
            {
                BaseLib::write_value_binary<unsigned int> (fin, sizeof(double)
                                                  * (long)  msh.getNumberOfElements());
                for (long j = 0; j < (long) msh.getNumberOfElements(); j++) {
                    double v = .0;
                    ele_values[i].second->eval(NumLib::TXPosition(j), v);
                    BaseLib::write_value_binary(fin, v);
                }
            }
        }
        else
            offset += (long) msh.getNumberOfElements() * sizeof(double)
                      + SIZE_OF_BLOCK_LENGTH_TAG;
        if (!_useBinary || !output_data)
            writeDataArrayFooter(fin);
    }

    return true;
}

unsigned char VtuWriter::getVTKCellType(const MeshLib::ElementShape::type ele_type)
{
    unsigned char cell_type = 0;

    switch(ele_type)
    {
    case MeshLib::ElementShape::LINE:               // vtk_line=3
        cell_type = 3;
        break;
    case MeshLib::ElementShape::QUAD:               // quadrilateral=9
        cell_type = 9;
        break;
    case MeshLib::ElementShape::HEXAHEDRON:         // hexahedron=12
        cell_type = 12;
        break;
    case MeshLib::ElementShape::TRIANGLE:           // triangle=5
        cell_type = 5;
        break;
    case MeshLib::ElementShape::TETRAHEDRON:        // tetrahedron=10
        cell_type = 10;
        break;
    case MeshLib::ElementShape::PRISM:              // wedge=13
        cell_type = 13;
        break;
    case MeshLib::ElementShape::PYRAMID:              // pyramid=14
        cell_type = 14;
        break;
    default:
        std::cerr << "***ERROR: NO CORRESPONDING VTK CELL TYPE FOUND. (ELEMENT TYPE=" <<
        ele_type << ")" << std::endl;
        break;
    }
    return cell_type;
}

void VtuWriter::initialize()
{
    //if (useBinary) {
    //======================================================================
    //# Set machine dependent stuff
    //Data type
    if (sizeof(unsigned char) == 1)
        type_UChar = UInt8;
    else if (sizeof(unsigned char) == 2)
        type_UChar = UInt16;
    if (sizeof(int) == 4)
        type_Int = Int32;
    else if (sizeof(int) == 8)
        type_Int = Int64;
    if (sizeof(unsigned int) == 4)
        type_UInt = UInt32;
    else if (sizeof(unsigned int) == 8)
        type_UInt = UInt64;
    if (sizeof(long) == 4)
        type_Long = Int32;
    else if (sizeof(long) == 8)
        type_Long = Int64;
    if (sizeof(double) == 4)
        type_Double = Float32;
    else if (sizeof(double) == 8)
        type_Double = Float64;
    //
    SIZE_OF_BLOCK_LENGTH_TAG = sizeof(unsigned int);
    //Endian(byte order)
    _isLittleEndian = BaseLib::IsLittleEndian();
    //}

    this->_isInitialized = true;
}

bool VtuWriter::writeDataArrayHeader(std::fstream &fin,
                                VTK_XML_DATA_TYPE data_type,
                                const std::string &str_name,
                                int nr_components,
                                const std::string &str_format,
                                long offset)
{
    std::string str_data_type;
    switch (data_type)
    {
    case Int8: str_data_type = "Int8";
        break;
    case UInt8: str_data_type = "UInt8";
        break;
    case Int16: str_data_type = "Int16";
        break;
    case UInt16: str_data_type = "UInt16";
        break;
    case Int32: str_data_type = "Int32";
        break;
    case UInt32: str_data_type = "UInt32";
        break;
    case Int64: str_data_type = "Int64";
        break;
    case UInt64: str_data_type = "UInt64";
        break;
    case Float32: str_data_type = "Float32";
        break;
    case Float64: str_data_type = "Float64";
        break;
    }
    fin << "        <DataArray type=\"" << str_data_type << "\"";
    if (str_name != "")
        fin << " Name=\"" << str_name << "\"";
    if (nr_components > 1)
        fin << " NumberOfComponents=\"" << nr_components << "\"";
    fin << " format=\"" << str_format << "\"";
    if (_useBinary)
        fin << " offset=\"" << offset << "\" /";
    fin << ">" << std::endl;

    return true;
}

bool VtuWriter::writeDataArrayFooter(std::fstream &fin)
{
    if (!_useBinary)
        fin << "        </DataArray>" << std::endl;

    return true;
}



bool VtuWriter::writePoints(std::fstream &fin, bool output_data, IMesh& msh, long &offset)
{
    const size_t n_msh_nodes = msh.getNumberOfNodes();
    if (output_data)
    {
        if (!_useBinary)
            for (size_t i = 0; i < n_msh_nodes; i++)
            {
                const GeoLib::Point* pnt = msh.getNodeCoordinatesRef(i);
                fin << "          " << (*pnt)[0] << " " << (*pnt)[1] << " " << (*pnt)[2] <<
                        std::endl;
            }
        else
        {
            BaseLib::write_value_binary<unsigned int>(fin, sizeof(double) * 3 * n_msh_nodes);
            for (size_t i = 0; i < n_msh_nodes; i++)
            {
                const GeoLib::Point* pnt = msh.getNodeCoordinatesRef(i);
                BaseLib::write_value_binary(fin, (*pnt)[0]);
                BaseLib::write_value_binary(fin, (*pnt)[1]);
                BaseLib::write_value_binary(fin, (*pnt)[2]);
            }
        }
    }
    else if (_useBinary)
        offset += n_msh_nodes * sizeof(double) * 3 + SIZE_OF_BLOCK_LENGTH_TAG;

    return true;
}

bool VtuWriter::writeCellConnectivity(std::fstream &fin,
                                        bool output_data,
                                        IMesh& msh,
                                        long &offset,
                                        long &sum_ele_components)
{
    if (output_data)
    {
        MeshLib::IElement* ele = NULL;
        if (!_useBinary)
            for (long i = 0; i < (long)msh.getNumberOfElements(); i++)
            {
                ele = msh.getElemenet(i);
                fin << "          ";
                for (size_t j = 0; j < ele->getNumberOfNodes(); j++)
                    fin << ele->getNodeID(j) << " ";
                fin << std::endl;
            }
        else
        {
            BaseLib::write_value_binary<unsigned int>(fin, sizeof(long) * sum_ele_components);
            for (long i = 0; i < (long)msh.getNumberOfElements(); i++)
            {
                ele = msh.getElemenet(i);
                for (size_t j = 0; j < ele->getNumberOfNodes(); j++)
                    BaseLib::write_value_binary<long>(fin, ele->getNodeID(j));
            }
        }
    }
    else if (_useBinary)
    {
        sum_ele_components = 0;
        MeshLib::IElement* ele = NULL;
        for (size_t i = 0; i < msh.getNumberOfElements(); i++) {
            ele = msh.getElemenet(i);
            sum_ele_components += ele->getNumberOfNodes();
        }
        offset += sum_ele_components * sizeof(long) + SIZE_OF_BLOCK_LENGTH_TAG;
    }

    return true;
}

bool VtuWriter::writeCellOffset(std::fstream &fin, bool output_data, IMesh& msh, long &offset)
{
    if (output_data)
    {
        MeshLib::IElement* ele = NULL;

        if (!_useBinary)
        {
            fin << "          ";
            long ele_offset = 0;
            for (long i = 0; i < (long)msh.getNumberOfElements(); i++)
            {
                ele = msh.getElemenet(i);
                ele_offset += ele->getNumberOfNodes();
                fin << ele_offset << " ";
            }
            fin << std::endl;
        }
        else
        {
            //OK411
            BaseLib::write_value_binary<unsigned int>(fin,
                                             sizeof(long) * (long)msh.getNumberOfElements());
            long ele_offset = 0;
            for (long i = 0; i < (long)msh.getNumberOfElements(); i++)
            {
                ele = msh.getElemenet(i);
                ele_offset += ele->getNumberOfNodes();
                BaseLib::write_value_binary(fin, ele_offset);
            }
        }
    }
    else if (_useBinary)
    {
        offset += (long)msh.getNumberOfElements() * sizeof(long) + SIZE_OF_BLOCK_LENGTH_TAG;
    }

    return true;
}

bool VtuWriter::writeCellType(std::fstream &fin, bool output_data, IMesh& msh, long &offset)
{
    if (output_data)
    {
        MeshLib::IElement* ele = NULL;
        if (!_useBinary)
        {
            fin << "          ";
            for(long i = 0; i < (long)msh.getNumberOfElements(); i++)
            {
                ele = msh.getElemenet(i);
                fin << (int)this->getVTKCellType(ele->getShapeType()) << " ";
            }
            fin << std::endl;
        }
        else
        {
            BaseLib::write_value_binary<unsigned int>(
                    fin,
                    sizeof(unsigned char) *
                    (long)msh.getNumberOfElements());
            for(long i = 0; i < (long)msh.getNumberOfElements(); i++)
            {
                ele = msh.getElemenet(i);
                BaseLib::write_value_binary(fin, this->getVTKCellType(ele->getShapeType()));
            }
        }
    }
    else if (_useBinary)
    {
        offset += (long)msh.getNumberOfElements() * sizeof(unsigned char) +
                  SIZE_OF_BLOCK_LENGTH_TAG;
    }

    return true;
}

