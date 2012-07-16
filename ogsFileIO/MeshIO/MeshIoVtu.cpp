
#include "MeshIoVtu.h"

#include <fstream>

#include "BaseLib/SystemTools.h"
#include "BaseLib/FileTools.h"

using namespace std;
using namespace MeshLib;

const std::string INDEX_STR = "  ";

MeshIoVtu::MeshIoVtu(bool binary_mode)
: useBinary(binary_mode)
{
	InitializeVTU();
}

bool MeshIoVtu::WriteXMLUnstructuredGrid(const std::string &vtkfile,
                                    MeshLib::IMesh &msh)
{
	//-------------------------------------------------------------------------
	//# Setup file stream
	//-------------------------------------------------------------------------
	std::fstream fin;
	if (useBinary)
		fin.open(vtkfile.data(), std::ios::out | std::ios::binary);
	else
		fin.open(vtkfile.data(), std::ios::out);

	if (!fin.good())
	{
		std::cout << "***Warning: Cannot open the output file, " << vtkfile << std::endl;
		return false;
	}

	if (!useBinary)
	{
		fin.setf(std::ios::scientific,std::ios::floatfield);
		fin.precision(12);
	}

	//-------------------------------------------------------------------------
	//# Output
	//-------------------------------------------------------------------------
	long offset = 0;

	string str_format;
	if (!useBinary)
		str_format = "ascii";
	else
		str_format = "appended";
	bool data_out = !useBinary;

	//# Header
	fin << "<?xml version=\"1.0\"?>" << std::endl;
	fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
	if (!useBinary || isLittleEndian)
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
	WriteDataArrayHeader(fin, type_Double, "", 3, str_format, offset);
	WriteMeshNodes(fin, data_out, msh, offset);
	WriteDataArrayFooter(fin);
	fin << "      </Points>" << endl;
	//....................................................................
	// Elements
	fin << "      <Cells>" << endl;
	//connectivity
	WriteDataArrayHeader(fin, type_Long, "connectivity", 0, str_format, offset);
	long sum_ele_components = 0;
	WriteMeshElementConnectivity(fin, data_out, msh, offset, sum_ele_components);
	WriteDataArrayFooter(fin);
	//offset
	WriteDataArrayHeader(fin, type_Long, "offsets", 0, str_format, offset);
	WriteMeshElementOffset(fin, data_out, msh, offset);
	WriteDataArrayFooter(fin);
	//type
	WriteDataArrayHeader(fin, type_UChar, "types", 0, str_format, offset);
	WriteMeshElementType(fin, data_out, msh, offset);
	WriteDataArrayFooter(fin);
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
	if (useBinary)
	{
		fin << "  <AppendedData encoding=\"raw\">" << endl;
		fin << "    _";

		//Node
		this->WriteMeshNodes(fin, true, msh, offset);
		//Element
		//conncectivity
		this->WriteMeshElementConnectivity(fin, true, msh, offset, sum_ele_components);
		//offset
		this->WriteMeshElementOffset(fin, true, msh, offset);
		//type
		this->WriteMeshElementType(fin, true, msh, offset);
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


unsigned char MeshIoVtu::GetVTKCellType(const MeshLib::ElementShape::type ele_type)
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
	}
	return cell_type;
}

void MeshIoVtu::InitializeVTU()
{
	//if (useBinary) {
	//======================================================================
	//# Set machine dependent stuff
	//Data type
	if (sizeof(unsigned char) == 1)
		type_UChar = MeshIoVtu::UInt8;
	else if (sizeof(unsigned char) == 2)
		type_UChar = MeshIoVtu::UInt16;
	if (sizeof(int) == 4)
		type_Int = MeshIoVtu::Int32;
	else if (sizeof(int) == 8)
		type_Int = MeshIoVtu::Int64;
	if (sizeof(unsigned int) == 4)
		type_UInt = MeshIoVtu::UInt32;
	else if (sizeof(unsigned int) == 8)
		type_UInt = MeshIoVtu::UInt64;
	if (sizeof(long) == 4)
		type_Long = MeshIoVtu::Int32;
	else if (sizeof(long) == 8)
		type_Long = MeshIoVtu::Int64;
	if (sizeof(double) == 4)
		type_Double = MeshIoVtu::Float32;
	else if (sizeof(double) == 8)
		type_Double = MeshIoVtu::Float64;
	//
	SIZE_OF_BLOCK_LENGTH_TAG = sizeof(unsigned int);
	//Endian(byte order)
	isLittleEndian = BaseLib::IsLittleEndian();
	//}

	this->isInitialized = true;
}

bool MeshIoVtu::WriteDataArrayHeader(std::fstream &fin,
                                VTK_XML_DATA_TYPE data_type,
                                const std::string &str_name,
                                int nr_components,
                                const std::string &str_format,
                                long offset)
{
	std::string str_data_type;
	switch (data_type)
	{
	case MeshIoVtu::Int8: str_data_type = "Int8";
		break;
	case MeshIoVtu::UInt8: str_data_type = "UInt8";
		break;
	case MeshIoVtu::Int16: str_data_type = "Int16";
		break;
	case MeshIoVtu::UInt16: str_data_type = "UInt16";
		break;
	case MeshIoVtu::Int32: str_data_type = "Int32";
		break;
	case MeshIoVtu::UInt32: str_data_type = "UInt32";
		break;
	case MeshIoVtu::Int64: str_data_type = "Int64";
		break;
	case MeshIoVtu::UInt64: str_data_type = "UInt64";
		break;
	case MeshIoVtu::Float32: str_data_type = "Float32";
		break;
	case MeshIoVtu::Float64: str_data_type = "Float64";
		break;
	}
	fin << "        <DataArray type=\"" << str_data_type << "\"";
	if (str_name != "")
		fin << " Name=\"" << str_name << "\"";
	if (nr_components > 1)
		fin << " NumberOfComponents=\"" << nr_components << "\"";
	fin << " format=\"" << str_format << "\"";
	if (useBinary)
		fin << " offset=\"" << offset << "\" /";
	fin << ">" << std::endl;

	return true;
}

bool MeshIoVtu::WriteDataArrayFooter(std::fstream &fin)
{
	if (!useBinary)
		fin << "        </DataArray>" << std::endl;

	return true;
}



bool MeshIoVtu::WriteMeshNodes(std::fstream &fin, bool output_data, IMesh& msh, long &offset)
{
	const size_t n_msh_nodes = msh.getNumberOfNodes();
	if (output_data)
	{
		if (!useBinary)
			for (size_t i = 0; i < n_msh_nodes; i++)
			{
				const GeoLib::Point* pnt = msh.getNodeCoordinatesRef(i);
				fin << "          " << (*pnt)[0] << " " << (*pnt)[1] << " " << (*pnt)[2] <<
				endl;
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
	else if (useBinary)
		offset += n_msh_nodes * sizeof(double) * 3 + SIZE_OF_BLOCK_LENGTH_TAG;

	return true;
}

bool MeshIoVtu::WriteMeshElementConnectivity(std::fstream &fin,
                                        bool output_data,
                                        IMesh& msh,
                                        long &offset,
                                        long &sum_ele_components)
{
	if (output_data)
	{
		MeshLib::IElement* ele = NULL;
		if (!useBinary)
			for (long i = 0; i < (long)msh.getNumberOfElements(); i++)
			{
				ele = msh.getElemenet(i);
				fin << "          ";
				for (size_t j = 0; j < ele->getNumberOfNodes(); j++)
					fin << ele->getNodeID(j) << " ";
				fin << endl;
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
	else if (useBinary)
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

bool MeshIoVtu::WriteMeshElementOffset(std::fstream &fin, bool output_data, IMesh& msh, long &offset)
{
	if (output_data)
	{
		MeshLib::IElement* ele = NULL;

		if (!useBinary)
		{
			fin << "          ";
			long ele_offset = 0;
			for (long i = 0; i < (long)msh.getNumberOfElements(); i++)
			{
				ele = msh.getElemenet(i);
				ele_offset += ele->getNumberOfNodes();
				fin << ele_offset << " ";
			}
			fin << endl;
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
	else if (useBinary)
	{
		offset += (long)msh.getNumberOfElements() * sizeof(long) + SIZE_OF_BLOCK_LENGTH_TAG;
	}

	return true;
}

bool MeshIoVtu::WriteMeshElementType(std::fstream &fin, bool output_data, IMesh& msh, long &offset)
{
	if (output_data)
	{
		MeshLib::IElement* ele = NULL;
		if (!useBinary)
		{
			fin << "          ";
			for(long i = 0; i < (long)msh.getNumberOfElements(); i++)
			{
				ele = msh.getElemenet(i);
				fin << (int)this->GetVTKCellType(ele->getShapeType()) << " ";
			}
			fin << endl;
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
				BaseLib::write_value_binary(fin, this->GetVTKCellType(ele->getShapeType()));
			}
		}
	}
	else if (useBinary)
	{
		offset += (long)msh.getNumberOfElements() * sizeof(unsigned char) +
		          SIZE_OF_BLOCK_LENGTH_TAG;
	}

	return true;
}


bool MeshIoVtu::WriteElementGroupID(std::fstream &fin,
                             bool output_data,
                             IMesh& msh,
                             long &offset)
{
	string str_format;
	if (!useBinary)
		str_format = "ascii";
	else
		str_format = "appended";


	if (!useBinary || !output_data)
		WriteDataArrayHeader(fin, this->type_Int, "MatGroup", 0, str_format, offset);
	if (output_data)
	{
		MeshLib::IElement* ele = NULL;
		if (!useBinary)
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

	if (!useBinary || !output_data)
		WriteDataArrayFooter(fin);


	return true;
}
