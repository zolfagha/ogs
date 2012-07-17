
#include "OutputIoVtu.h"

#include <fstream>

#include "BaseLib/FileTools.h"
#include "MeshLib/Core/CoordinateSystem.h"
#include "MeshLib/Core/MeshGeometricProperties.h"

using namespace MeshLib;

const std::string INDEX_STR = "  ";
const std::string velocity_name[3][4] =
{
	{
		"VELOCITY_X1", "VELOCITY_Y1", "VELOCITY_Z1", "NODAL_VELOCITY1"
	}
	,
	{
		"VELOCITY_X2", "VELOCITY_Y2", "VELOCITY_Z2", "NODAL_VELOCITY2"
	}
	,
	{
		"VELOCITY1_X", "VELOCITY1_Y", "VELOCITY1_Z", "GL_NODAL_VELOCITY1"
	}
};

OutputIoVtu::OutputIoVtu(bool binary_mode)
: MeshIoVtu(binary_mode)
{
}

bool OutputIoVtu::Write(const std::string &vtkfile,
                                    MeshLib::IMesh &msh, std::vector<NodalValue> &nodal_values, std::vector<ElementalValue> &elemental_values)
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

	std::string str_format;
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
	fin << "  <UnstructuredGrid>" << std::endl;
	fin << "    <Piece NumberOfPoints=\"" << msh.getNumberOfNodes() <<
	"\" NumberOfCells=\"" << msh.getNumberOfElements() << "\">" << std::endl;
	//....................................................................
	// Nodes
	fin << "      <Points>" << std::endl;
	WriteDataArrayHeader(fin, type_Double, "", 3, str_format, offset);
	WriteMeshNodes(fin, data_out, msh, offset);
	WriteDataArrayFooter(fin);
	fin << "      </Points>" << std::endl;
	//....................................................................
	// Elements
	fin << "      <Cells>" << std::endl;
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
	fin << "      </Cells>" << std::endl;
	//....................................................................
	// Nodal values
	if (nodal_values.size() > 0)
		fin << "      <PointData Scalars=\"" << nodal_values[0].first << "\">" << std::endl;
	else
		fin << "      <PointData Scalars=\"scalars\">" << std::endl;
	WriteNodalValue(fin, data_out, nodal_values, msh, offset);
	fin << "      </PointData>" << std::endl;
	//....................................................................
	// Element values
	fin << "      <CellData>" << std::endl;
	WriteElementValue(fin, data_out, elemental_values, msh, offset);
	fin << "      </CellData>" << std::endl;
	fin << "    </Piece>" << std::endl;
	fin << "  </UnstructuredGrid>" << std::endl;

	//======================================================================
	// Raw data (for binary mode)
	if (useBinary)
	{
		fin << "  <AppendedData encoding=\"raw\">" << std::endl;
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
		this->WriteNodalValue(fin, true, nodal_values, msh, offset);
		// Elemental values
		this->WriteElementValue(fin, true, elemental_values, msh, offset);

		fin << std::endl;
		fin << "  </AppendedData>" << std::endl;
	}

	fin << "</VTKFile>" << std::endl;
	fin.close();

	return true;
}

bool OutputIoVtu::WriteNodalValue(std::fstream &fin,
                           bool output_data,
                           std::vector<NodalValue> &nod_values,
                           IMesh &msh,
                           long &offset)
{
	std::vector<int> NodeIndex(nod_values.size());

	std::string str_format;
	if (!this->useBinary)
		str_format = "ascii";
	else
		str_format = "appended";

    bool isXZplane = (msh.getGeometricProperty()->getCoordinateSystem()->getType() == MeshLib::CoordinateSystemType::XZ);
    bool is3D = (msh.getDimension() == 3);
	bool outNodeVelocity = false;
    bool outNodeDisplacement = false;

	//Nodal values
	for (int i = 0; i < (int) nod_values.size(); i++)
	{
        const std::string &var_name = nod_values[i].first;
		//is velocity
		if (var_name.find("VELOCITY") != std::string::npos)
		{
			outNodeVelocity = true;
			continue;
		}
        if (var_name.find("DISPLACEMENT") != std::string::npos)
        {
            outNodeDisplacement = true;
            continue;
        }

		if (!useBinary || !output_data)
			WriteDataArrayHeader(fin, type_Double, var_name, 0, str_format, offset);

		if (output_data)
		{
			if (!useBinary) {
				fin << "          ";
			} else {
				BaseLib::write_value_binary<unsigned int> (fin, sizeof(double)* msh.getNumberOfNodes());
			}
            for (size_t j = 0; j < msh.getNumberOfNodes(); j++) {
                double v = nod_values[i].second[j];
                if (!useBinary) {
                    fin << v << " ";
                } else {
                    BaseLib::write_value_binary(fin, v);
                }
            }
            if (!useBinary) {
                fin << std::endl;
            }
		}
		else
			offset += msh.getNumberOfNodes() * sizeof(double)
			          + SIZE_OF_BLOCK_LENGTH_TAG;

		if (!useBinary || !output_data)
			WriteDataArrayFooter(fin);
	}

#if 0
	// Nodal velocities
	if (outNodeVelocity)
	{
		unsigned int velocity_id = 0;
		for (int i = 0; i < (int) nod_values.size(); i++)
		{
            const string &internal_val_name = nod_values[i].first;
//            const string &external_val_name = out->_alias_nod_value_vector[i];
			if (internal_val_name.find("VELOCITY_X1") != string::npos)
			{
				velocity_id = 0;
			}
			else if (internal_val_name.find("VELOCITY_X2")
			         != string::npos)
			{
				velocity_id = 1;
			}
			else if (internal_val_name.find("VELOCITY1_X")
			         != string::npos)
			{
				velocity_id = 2;
			}
			else
				continue;

			if (!useBinary || !output_data)
				WriteDataArrayHeader(fin,
				                     this->type_Double,
				                     velocity_name[velocity_id][3],
				                     3,
				                     str_format,
				                     offset);
			if (output_data)
			{
				int ix, iy, iz;
				ix = nod_values[i].second(velocity_name[velocity_id][0],true); // JT: Fix. Need latest value.
				iy = m_pcs->GetNodeValueIndex(velocity_name[velocity_id][1],true);
				iz = m_pcs->GetNodeValueIndex(velocity_name[velocity_id][2],true);
				if (!useBinary)
				{
					fin << "          ";
					for (size_t j = 0l; j < msh.getNumberOfNodes(); j++)
					{
						fin << m_pcs->GetNodeValue(
						        msh->nod_vector[j]->GetIndex(), ix) << " ";
                        if (!isXZplane) {
                            fin << m_pcs->GetNodeValue(
                                msh->nod_vector[j]->GetIndex(), iy) << " ";
                            fin << m_pcs->GetNodeValue(
                                msh->nod_vector[j]->GetIndex(), iz) << " ";
                        } else {
                            fin << m_pcs->GetNodeValue(
                                msh->nod_vector[j]->GetIndex(), iz) << " ";
                            fin << m_pcs->GetNodeValue(
                                msh->nod_vector[j]->GetIndex(), iy) << " ";
                        }
					}
					fin << endl;
				}
				else
				{
					BaseLib::write_value_binary<unsigned int> (fin, sizeof(double)
					                                  * msh->GetNodesNumber(
					                                          false) * 3);
					for (size_t j = 0l; j < msh.getNumberOfNodes(); j++)
					{
						BaseLib::write_value_binary(fin, m_pcs->GetNodeValue(
						                           msh->nod_vector[j]->
						                           GetIndex(), ix));
                        if (!isXZplane) {
                            BaseLib::write_value_binary(fin, m_pcs->GetNodeValue(
                                msh->nod_vector[j]->
                                GetIndex(), iy));
                            BaseLib::write_value_binary(fin, m_pcs->GetNodeValue(
                                msh->nod_vector[j]->
                                GetIndex(), iz));
                        } else {
                            BaseLib::write_value_binary(fin, m_pcs->GetNodeValue(
                                msh->nod_vector[j]->
                                GetIndex(), iz));
                            BaseLib::write_value_binary(fin, m_pcs->GetNodeValue(
                                msh->nod_vector[j]->
                                GetIndex(), iy));
                        }
					}
				}
			}
			else {
				offset += msh.getNumberOfNodes() * 3 * sizeof( double) +
				          SIZE_OF_BLOCK_LENGTH_TAG;
			}

			if (!useBinary || !output_data)
				WriteDataArrayFooter(fin);
		}
	}
#endif


    return true;
}

bool OutputIoVtu::WriteElementValue(std::fstream &fin,
                             bool output_data,
                             std::vector<ElementalValue> &ele_values,
                             IMesh &msh,
                             long &offset)
{
	std::vector<int> ele_value_index_vector(ele_values.size());

    bool isXZplane = (msh.getGeometricProperty()->getCoordinateSystem()->getType() == MeshLib::CoordinateSystemType::XZ);

    std::string str_format;
	if (!this->useBinary)
		str_format = "ascii";
	else
		str_format = "appended";

	bool outEleVelocity = false;

	//Element values
	for (int i = 0; i < (int) ele_value_index_vector.size(); i++)
	{
		if (ele_values[i].first.find("VELOCITY") != std::string::npos)
		{
			outEleVelocity = true;
			continue;
		}

		if (!useBinary || !output_data)
			WriteDataArrayHeader(fin, this->type_Double,
					ele_values[i].first, 0, str_format, offset);
		if (output_data)
		{
			if (!useBinary)
			{
				fin << "          ";
				for (long j = 0; j < (long) msh.getNumberOfElements(); j++)
					fin << ele_values[i].second[j] << " ";
				fin << std::endl;
			}
			else
			{
				BaseLib::write_value_binary<unsigned int> (fin, sizeof(double)
				                                  * (long)  msh.getNumberOfElements());
				for (long j = 0; j < (long) msh.getNumberOfElements(); j++)
					BaseLib::write_value_binary(fin, ele_values[i].second[j]);
			}
		}
		else
			offset += (long) msh.getNumberOfElements() * sizeof(double)
			          + SIZE_OF_BLOCK_LENGTH_TAG;  //OK411
		if (!useBinary || !output_data)
			WriteDataArrayFooter(fin);
	}

	return true;
}

