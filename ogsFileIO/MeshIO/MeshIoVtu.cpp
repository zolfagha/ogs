
#include "MeshIoVtu.h"

#include <fstream>

#include "BaseLib/SystemTools.h"
#include "BaseLib/FileTools.h"

using namespace std;
using namespace MeshLib;

MeshIoVtu::MeshIoVtu(bool binary_mode)
: VtuWriter(binary_mode)
{
}

bool MeshIoVtu::WriteXMLUnstructuredGrid(const std::string &vtkfile,
                                    MeshLib::IMesh &msh)
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

    string str_format;
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
    fin << "  <UnstructuredGrid>" << endl;
    fin << "    <Piece NumberOfPoints=\"" << msh.getNumberOfNodes() <<
    "\" NumberOfCells=\"" << msh.getNumberOfElements() << "\">" << std::endl;
    //....................................................................
    // Nodes
    fin << "      <Points>" << std::endl;
    writeDataArrayHeader(fin, type_Double, "", 3, str_format, offset);
    writePoints(fin, data_out, msh, offset);
    writeDataArrayFooter(fin);
    fin << "      </Points>" << endl;
    //....................................................................
    // Elements
    fin << "      <Cells>" << endl;
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
    fin << "      </Cells>" << endl;
    //....................................................................
    // Element values
    fin << "      <CellData>" << endl;
    WriteElementGroupID(fin, data_out, msh, offset);
    fin << "      </CellData>" << endl;
    fin << "    </Piece>" << endl;
    fin << "  </UnstructuredGrid>" << endl;

    //======================================================================
    // Raw data (for binary mode)
    if (_useBinary)
    {
        fin << "  <AppendedData encoding=\"raw\">" << endl;
        fin << "    _";

        //Node
        this->writePoints(fin, true, msh, offset);
        //Element
        //conncectivity
        this->writeCellConnectivity(fin, true, msh, offset, sum_ele_components);
        //offset
        this->writeCellOffset(fin, true, msh, offset);
        //type
        this->writeCellType(fin, true, msh, offset);
        // Nodal values
        //this->WriteNodalValue(fin, true, out, msh, offset);
        // Elemental values
        this->WriteElementGroupID(fin, true,  msh, offset);

        fin << endl;
        fin << "  </AppendedData>" << endl;
    }

    fin << "</VTKFile>" << endl;
    fin.close();

    return true;
}



bool MeshIoVtu::WriteElementGroupID(std::fstream &fin,
                             bool output_data,
                             IMesh& msh,
                             long &offset)
{
    string str_format;
    if (!_useBinary)
        str_format = "ascii";
    else
        str_format = "appended";


    if (!_useBinary || !output_data)
        writeDataArrayHeader(fin, this->type_Int, "MatGroup", 0, str_format, offset);
    if (output_data)
    {
        MeshLib::IElement* ele = NULL;
        if (!_useBinary)
        {
            fin << "          ";
            for(long i = 0; i < (long)msh.getNumberOfElements(); i++)
            {
                ele = msh.getElemenet(i);
                fin << ele->getGroupID() << " ";
            }
            fin << endl;
        }
        else
        {
            BaseLib::write_value_binary<unsigned int>(
                    fin,
                    sizeof(int) *
                    (long)msh.getNumberOfElements());
            for (long i = 0; i < (long)msh.getNumberOfElements(); i++)
                BaseLib::write_value_binary(fin, msh.getElemenet(i)->getGroupID());
        }
    }
    else
    {
        offset += (long)msh.getNumberOfElements() * sizeof(int) +
                  SIZE_OF_BLOCK_LENGTH_TAG;
    }

    if (!_useBinary || !output_data)
        writeDataArrayFooter(fin);


    return true;
}
