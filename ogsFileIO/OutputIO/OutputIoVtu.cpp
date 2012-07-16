
#include "OutputIoVtu.h"

#include <fstream>

#include "BaseLib/FileTools.h"

using namespace std;
using namespace MeshLib;

const std::string INDEX_STR = "  ";

OutputIoVtu::OutputIoVtu(bool binary_mode)
: MeshIoVtu(binary_mode)
{
}

bool OutputIoVtu::WriteXMLUnstructuredGrid(const std::string &vtkfile,
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

bool OutputIoVtu::WriteNodalValue(std::fstream &fin,
                           bool output_data,
                           std::vector<NodalValue> &nod_values,
                           IMesh &msh,
                           long &offset)
{
	std::vector<int> NodeIndex(nod_values.size());

 	string str_format;
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
        const string &var_name = nod_values[i].first;
		//is velocity
		if (var_name.find("VELOCITY") != string::npos)
		{
			outNodeVelocity = true;
			continue;
		}
        if (var_name.find("DISPLACEMENT") != string::npos)
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
                fin << endl;
            }
		}
		else
			offset += msh.getNumberOfNodes() * sizeof(double)
			          + SIZE_OF_BLOCK_LENGTH_TAG;

		if (!useBinary || !output_data)
			WriteDataArrayFooter(fin);
	}

	// Nodal velocities
	if (outNodeVelocity)
	{
		unsigned int velocity_id = 0;
		for (int i = 0; i < (int) out->_nod_value_vector.size(); i++)
		{
            const string &internal_val_name = out->_nod_value_vector[i];
//            const string &external_val_name = out->_alias_nod_value_vector[i];
			if (internal_val_name.find("VELOCITY_X1") != string::npos)
			{
				if (out->m_pcs == NULL)
					m_pcs = PCSGet(internal_val_name, true);
				velocity_id = 0;
			}
			else if (internal_val_name.find("VELOCITY_X2")
			         != string::npos)
			{
				if (out->m_pcs == NULL)
					m_pcs = PCSGet(internal_val_name, true);
				velocity_id = 1;
			}
			else if (internal_val_name.find("VELOCITY1_X")
			         != string::npos)
			{
				if (out->m_pcs == NULL)
					m_pcs = PCSGet(internal_val_name, true);
				velocity_id = 2;
			}
			else
				continue;
			if (!m_pcs)
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
				ix = m_pcs->GetNodeValueIndex(velocity_name[velocity_id][0],true); // JT: Fix. Need latest value.
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
			else
				offset += msh.getNumberOfNodes() * 3 * sizeof( double) +
				          SIZE_OF_BLOCK_LENGTH_TAG;

			if (!useBinary || !output_data)
				WriteDataArrayFooter(fin);
		}
	}

    //Displacement
    if (outNodeDisplacement)
    {
//        unsigned int disp_id = 0;
        for (int i = 0; i < (int)out->_nod_value_vector.size(); i++)
        {
            const string &internal_val_name = out->_nod_value_vector[i];
//            const string &external_val_name = out->_alias_nod_value_vector[i];
            if (internal_val_name.find("DISPLACEMENT_X1") != string::npos)
            {
                if (out->m_pcs == NULL)
                    m_pcs = PCSGet(internal_val_name,true);
//                disp_id = 0;
            }
            else
                continue;
            if(!m_pcs)
                continue;

            if (!useBinary || !output_data)
                WriteDataArrayHeader(fin, this->type_Double, "DISPLACEMENT", 3, str_format, offset);
            if (output_data)
            {
                int var_id[3] = {};
                var_id[0] = m_pcs->GetNodeValueIndex("DISPLACEMENT_X1");
                var_id[1] = m_pcs->GetNodeValueIndex("DISPLACEMENT_Y1");
                var_id[2] = -1;
                if (is3D) {
                    var_id[2] = m_pcs->GetNodeValueIndex("DISPLACEMENT_Z1");
                } else if (isXZplane) {
                    var_id[1] = -1;
                    var_id[2] = m_pcs->GetNodeValueIndex("DISPLACEMENT_Y1");
                }
                //
                if (!useBinary) {
                    fin << "          ";
                } else {
                	BaseLib::write_value_binary<unsigned int>(fin, sizeof(double)*msh.getNumberOfNodes()*3);
                }
                double u[3] = {};
                for(size_t j = 0l; j < msh.getNumberOfNodes(); j++)
                {
                    for (size_t k=0; k<3; k++) {
                        if (var_id[k]<0)
                            u[k] = .0;
                        else
                            u[k] =  m_pcs->GetNodeValue(msh->nod_vector[j]->GetIndex(), var_id[k]);
                    }

                    if (!useBinary) {
                        for (size_t k=0; k<3; k++)
                            fin << u[k] << " ";
                    } else {
                        for (size_t k=0; k<3; k++)
                            BaseLib::write_value_binary(fin, u[k]);
                    }
                }
                if (!useBinary) {
                    fin << endl;
                } else {
                    BaseLib::write_value_binary<unsigned int>(fin, sizeof(double)*msh.getNumberOfNodes()*3);
                }
            }
            else
                offset += msh.getNumberOfNodes() * 3 * sizeof(double) +
                SIZE_OF_BLOCK_LENGTH_TAG;
            if (!useBinary || !output_data)
                WriteDataArrayFooter(fin);
        }
    }

    return true;
}

bool OutputIoVtu::WriteElementValue(std::fstream &fin,
                             bool output_data,
                             std::vector<ElementalValue> &ele_values,
                             IMesh &msh,
                             long &offset)
{
	std::vector<int> ele_value_index_vector(out->getElementValueVector().size());
	if (ele_value_index_vector.size() > 0) // GetELEValuesIndexVector() should check this!
		out->GetELEValuesIndexVector(ele_value_index_vector);
	CRFProcess* m_pcs = NULL;
	MeshLib::CElem* ele = NULL;

    bool isXZplane = (msh->GetCoordinateFlag()==22);

	string str_format;
	if (!this->useBinary)
		str_format = "ascii";
	else
		str_format = "appended";

	bool outEleVelocity = false;

	//Element values
	for (int i = 0; i < (int) ele_value_index_vector.size(); i++)
	{
		if (out->getElementValueVector()[i].find("VELOCITY") != string::npos)
		{
			outEleVelocity = true;
			continue;
		}
		m_pcs = out->GetPCS_ELE(out->getElementValueVector()[i]);

		if (!useBinary || !output_data)
			WriteDataArrayHeader(fin, this->type_Double,
			                     out->getElementValueVector()[i], 0, str_format, offset);
		if (output_data)
		{
			if (!useBinary)
			{
				fin << "          ";
				for (long j = 0; j < (long) msh->ele_vector.size(); j++)
					fin << m_pcs->GetElementValue(j, ele_value_index_vector[i])
					    << " ";
				fin << endl;
			}
			else
			{
				BaseLib::write_value_binary<unsigned int> (fin, sizeof(double)
				                                  * (long) msh->ele_vector.size()); //OK411
				for (long j = 0; j < (long) msh->ele_vector.size(); j++)
					BaseLib::write_value_binary(fin, m_pcs->GetElementValue(j,
					                                               ele_value_index_vector
					                                               [i]));
			}
		}
		else
			offset += (long) msh->ele_vector.size() * sizeof(double)
			          + SIZE_OF_BLOCK_LENGTH_TAG;  //OK411
		if (!useBinary || !output_data)
			WriteDataArrayFooter(fin);
	}

	//Element velocity
	if (outEleVelocity)
	{
		if (!useBinary || !output_data)
			WriteDataArrayHeader(fin, this->type_Double, "ELEMENT_VELOCITY", 3,
			                     str_format, offset);
		if (output_data)
		{
			if (!useBinary)
			{
				fin << "          ";
				static double ele_vel[3] = { 0.0, 0.0, 0.0 };
				for (long i = 0; i < (long) msh->ele_vector.size(); i++)
				{
					ele_gp_value[i]->getIPvalue_vec(0, ele_vel);
					fin << ele_vel[0] << " ";
                    if (!isXZplane) {
                        fin << ele_vel[1] << " ";
                        fin << ele_vel[2] << " ";
                    } else {
                        fin << ele_vel[2] << " ";
                        fin << ele_vel[1] << " ";
                    }
				}
				fin << endl;
			}
			else
			{
				static double ele_vel[3] = { 0.0, 0.0, 0.0 };
				BaseLib::write_value_binary<unsigned int> (fin, sizeof(double) * 3 *
				                                  (long )msh->ele_vector.size()); //OK411
				for(long i = 0; i < (long)msh->ele_vector.size(); i++)
				{
					ele_gp_value[i]->getIPvalue_vec(0, ele_vel);
                    BaseLib::write_value_binary(fin, ele_vel[0]);
                    if (!isXZplane) {
                        BaseLib::write_value_binary(fin, ele_vel[1]);
                        BaseLib::write_value_binary(fin, ele_vel[2]);
                    } else {
                        BaseLib::write_value_binary(fin, ele_vel[2]);
                        BaseLib::write_value_binary(fin, ele_vel[1]);
                    }
				}
			}
		}
		else
			//OK411
			offset += (long)msh->ele_vector.size() * sizeof(double) * 3 +
			          SIZE_OF_BLOCK_LENGTH_TAG;
		if (!useBinary || !output_data)
			WriteDataArrayFooter(fin);

		//    if(out->pcs_type_name.compare("FLUID_MOMENTUM")==0)
		if(out->getProcessType () == FiniteElement::FLUID_MOMENTUM)
		{
			if (!useBinary || !output_data)
				WriteDataArrayHeader(fin,
				                     this->type_Double,
				                     "GLOBAL_VELOCITY",
				                     3,
				                     str_format,
				                     offset);
			if (output_data)
			{
				CRFProcess* pch_pcs = PCSGet("FLUID_MOMENTUM");
				if (!this->useBinary)
				{
					fin << "          ";
					for(long i = 0; i < (long)msh->ele_vector.size(); i++)
					{
						fin << pch_pcs->GetElementValue(
						        i,
						        pch_pcs->
						        GetElementValueIndex("VELOCITY1_X") +
						        1) << " ";
                        if (!isXZplane) {
                            fin << pch_pcs->GetElementValue(
                                i,
                                pch_pcs->
                                GetElementValueIndex("VELOCITY1_Y") +
                                1) << " ";
                            fin << pch_pcs->GetElementValue(
                                i,
                                pch_pcs->
                                GetElementValueIndex("VELOCITY1_Z") +
                                1) << " ";
                        } else {
                            fin << pch_pcs->GetElementValue(
                                i,
                                pch_pcs->
                                GetElementValueIndex("VELOCITY1_Z") +
                                1) << " ";
                            fin << pch_pcs->GetElementValue(
                                i,
                                pch_pcs->
                                GetElementValueIndex("VELOCITY1_Y") +
                                1) << " ";
                        }
					}
					fin << endl;
				}
				else
				{
					//OK411
					BaseLib::write_value_binary<unsigned int>(fin, sizeof(double) * 3 *
					                                 (long)msh->ele_vector.size());
					for(long i = 0; i < (long)msh->ele_vector.size(); i++)
					{
						BaseLib::write_value_binary(fin,
						                   pch_pcs->GetElementValue(i,
						                                            pch_pcs
						                                            ->
						                                            GetElementValueIndex(
						                                                    "VELOCITY1_X")
						                                            + 1));
                        if (!isXZplane) {
                            BaseLib::write_value_binary(fin,
                                pch_pcs->GetElementValue(i,
                                pch_pcs
                                ->
                                GetElementValueIndex(
                                "VELOCITY1_Y")
                                + 1));
                            BaseLib::write_value_binary(fin,
                                pch_pcs->GetElementValue(i,
                                pch_pcs
                                ->
                                GetElementValueIndex(
                                "VELOCITY1_Z")
                                + 1));
                        } else {
                            BaseLib::write_value_binary(fin,
                                pch_pcs->GetElementValue(i,
                                pch_pcs
                                ->
                                GetElementValueIndex(
                                "VELOCITY1_Z")
                                + 1));
                        }
					}
				}
			}
			else
				//OK411
				offset += (long)msh->ele_vector.size() * sizeof(double) * 3 +
				          SIZE_OF_BLOCK_LENGTH_TAG;
			if (!useBinary || !output_data)
				WriteDataArrayFooter(fin);
		}
	}
	//Material information
	if(mmp_vector.size() > 1)
	{
		if (!useBinary || !output_data)
			WriteDataArrayHeader(fin, this->type_Int, "MatGroup", 0, str_format, offset);
		if (output_data)
		{
			if (!this->useBinary)
			{
				fin << "          ";
				for(long i = 0; i < (long)msh->ele_vector.size(); i++)
				{
					ele = msh->ele_vector[i];
					fin << ele->GetPatchIndex() << " ";
				}
				fin << endl;
			}
			else
			{
				//OK411
				BaseLib::write_value_binary<unsigned int>(
				        fin,
				        sizeof(int) *
				        (long)msh->ele_vector.size());
				for (long i = 0; i < (long)msh->ele_vector.size(); i++)
					BaseLib::write_value_binary(fin, msh->ele_vector[i]->GetPatchIndex());
			}
		}
		else
			//OK411
			offset += (long)msh->ele_vector.size() * sizeof(int) +
			          SIZE_OF_BLOCK_LENGTH_TAG;
		if (!useBinary || !output_data)
			WriteDataArrayFooter(fin);
	}
	return true;
}

