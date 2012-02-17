/**************************************************************************
   PARLib - Object:
   Task:
   Programing:
   07/2004 OK Implementation
   07/2004 OK Version 1 untill 3.9.17OK6
   07/2004 OK Version 2 from   3.9.17OK7
   07/2006 WW Local topology, High order nodes
   last modified:
**************************************************************************/
//---- MPI Parallel --------------
// There is a name conflict between stdio.h and the MPI C++ binding
// with respect to the names SEEK_SET, SEEK_CUR, and SEEK_END.  MPI
// wants these in the MPI namespace, but stdio.h will #define these
// to integer values.  #undef'ing these can cause obscure problems
// with other include files (such as iostream), so we instead use
// #error to indicate a fatal error.  Users can either #undef
// the names before including mpi.h or include mpi.h *before* stdio.h
// or iostream.
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || \
        defined(USE_MPI_GEMS)
//#undef SEEK_SET  //WW
//#undef SEEK_END  //WW
//#undef SEEK_CUR  //WW
#include <mpi.h>
int size;
int myrank;
int mysize;
char t_fname[3];
double time_ele_paral;
#endif
//---- MPI Parallel --------------

#include <cmath>
#include <iostream>

#include "Base/CodingTools.h"

#include "par_ddc.h"
#include "rf_num_new.h"
#include "equation_class.h"
#include "matrix_class.h"

#define MAX_ZEILE 512
using namespace std;


namespace NumLib
{
namespace OGS5
{


bool KeywordFound(const string &line)
{
	string hash("#");
	if(line.find(hash) != string::npos)
		return true;
	else
		return false;
}

bool SubKeywordFound(const string &line)
{
	string dollar("$");
	if(line.find(dollar) != string::npos)
		return true;
	else
		return false;
}





//////////////////////////////////////////////////////////////////////////
// CPARDomain

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   07/2004 OK Implementation
**************************************************************************/
CPARDomain::CPARDomain(void)
{
#if defined(USE_MPI)                           // 13.12.2007 WW
	t_border_nodes = NULL;
	t_border_nodes_size = t_border_nodes_sizeH = 0;
	//
#if defined(NEW_BREDUCE)
	receive_cnt_b = new int[mysize];
	receive_disp_b = new int[mysize];
#endif
	receive_cnt_i = new int[mysize];
	receive_disp_i = new int[mysize];
	receive_cnt = new int[mysize];
	receive_disp = new int[mysize];
#endif
}

CPARDomain::~CPARDomain(void)
{
	elements.clear();
	nodes.clear();
	nodes_inner.clear();
	nodes_halo.clear();
	long i;

	// WW
	for(i = 0; i < (long)element_nodes_dom.size(); i++)
	{
		delete [] element_nodes_dom[i];
		element_nodes_dom[i] = NULL;
	}
	for(long i = 0; i < (long)node_conneted_nodes.size(); i++)
	{
		delete [] node_conneted_nodes[i];
		node_conneted_nodes[i] = NULL;
	}

	//
    Base::releaseObjectsInStdVector(_vec_sparse);
    Base::releaseObjectsInStdVector(_vec_eqs);

#if defined(USE_MPI)                           // 13.12.2007 WW
	//
	if(t_border_nodes)
		delete [] t_border_nodes;
	t_border_nodes = NULL;
#if defined(NEW_BREDUCE)
	delete [] receive_cnt_b;
	delete [] receive_disp_b;
	receive_cnt_b = NULL;
	receive_disp_b = NULL;
#endif
	delete [] receive_cnt_i;
	delete [] receive_disp_i;
	delete [] receive_cnt;
	delete [] receive_disp;
	receive_cnt_i = NULL;
	receive_disp_i = NULL;
	receive_cnt = NULL;
	receive_disp = NULL;
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task: Release partial memory
   Programing:
   012/2007 WW Implementation
**************************************************************************/
#if defined(USE_MPI)
void CPARDomain::ReleaseMemory()
{
	// WW
	for(long i = 0; i < (long)element_nodes_dom.size(); i++)
	{
		delete [] element_nodes_dom[i];
		element_nodes_dom[i] = NULL;
	}
	for(long i = 0; i < (long)node_conneted_nodes.size(); i++)
	{
		delete [] node_conneted_nodes[i];
		node_conneted_nodes[i] = NULL;
	}
	if(eqs)
		delete eqs;
	if(eqsH)
		delete eqsH;
	if(sparse_graph)
		delete sparse_graph;
	if(sparse_graph_H)
		delete sparse_graph_H;
	if(t_border_nodes)
		delete [] t_border_nodes;
	t_border_nodes = NULL;
#if defined(NEW_BREDUCE)
	delete [] receive_cnt_b;
	delete [] receive_disp_b;
	receive_cnt_b = NULL;
	receive_disp_b = NULL;
#endif
	delete [] receive_cnt_i;
	delete [] receive_disp_i;
	delete [] receive_cnt;
	delete [] receive_disp;
	receive_cnt_i = NULL;
	receive_disp_i = NULL;
	receive_cnt = NULL;
	receive_disp = NULL;
}
#endif
/**************************************************************************
   FEMLib-Method:
   Task: ST read function
   Programing:
   07/2004 OK Implementation
   07/2004 OK Version 1 untill 3.9.17OK6
   07/2004 OK Version 2 from   3.9.17OK7
**************************************************************************/
ios::pos_type CPARDomain::Read(ifstream* ddc_file)
{
	char line[MAX_ZEILE];
	string sub_line;
	string sub_string;
	string cut_string;
	string line_string;
	string delimiter(" ");
	string delimiter_type(";");
	bool new_subkeyword = false;
	string dollar("$");
	bool new_keyword = false;
	string hash("#");
	ios::pos_type position;
	long i;
	//  CElem* m_ele = NULL;
	//======================================================================
	while (!new_keyword)
	{
		position = ddc_file->tellg();
		ddc_file->getline(line,MAX_ZEILE);
		line_string = line;
		if(line_string.find(hash) != string::npos)
		{
			new_keyword = true;
			break;
		}
		//....................................................................
		// subkeyword found
		if(line_string.find("$ELEMENTS") != string::npos)
			while ((!new_keyword) && (!new_subkeyword))
			{
				position = ddc_file->tellg();
				ddc_file->getline(line,MAX_ZEILE);
				line_string = line;
				if(line_string.find(hash) != string::npos)
				{
					new_keyword = true;
					break;
				}
				if(line_string.find(dollar) != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				i = strtol(line,NULL,0);
				elements.push_back(i);
			}
		//....................................................................
		// subkeyword found
		if(line_string.find("$NODES_INNER") != string::npos)
			while (!new_keyword)
			{
				position = ddc_file->tellg();
				ddc_file->getline(line,MAX_ZEILE);
				line_string = line;
				if(line_string.find(hash) != string::npos)
				{
					new_keyword = true;
					break;
				}
				if(line_string.find(dollar) != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				i = strtol(line,NULL,0);
				nodes_inner.push_back(i);
			}
		//....................................................................
		// subkeyword found
		if(line_string.find("$NODES_BORDER") != string::npos)
			while (!new_keyword)
			{
				position = ddc_file->tellg();
				ddc_file->getline(line,MAX_ZEILE);
				line_string = line;
				if(line_string.find(hash) != string::npos)
				{
					new_keyword = true;
					break;
				}
				if(line_string.find(dollar) != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				i = strtol(line,NULL,0);
				nodes_halo.push_back(i);
			}
		//....................................................................
	}
	//======================================================================
	return position;
}

/**************************************************************************
   PARLib-Method:
   Task:
   Programing:
   07/2004 OK Implementation
**************************************************************************/
void CPARDomain::CreateNodes()
{
	long j,k;
	long no_nodes_halo, no_nodes_inner;
	//----------------------------------------------------------------------
	no_nodes_inner = (long)nodes_inner.size();
	for(j = 0; j < no_nodes_inner; j++)
	{
		k = nodes_inner[j];
		nodes.push_back(k);
		//cout << nodes[j] << endl;
	}
	//cout << "---" << endl;
	no_nodes_halo = (long)nodes_halo.size();
	for(j = 0; j < no_nodes_halo; j++)
	{
		k = nodes_halo[j];
		nodes.push_back(k);
		//cout << nodes[no_nodes_inner+j] << endl;
	}
	nnodes_dom = no_nodes_halo + no_nodes_inner; //WW
	nnodesHQ_dom = nnodes_dom;

	//----------------------------------------------------------------------
}

/**************************************************************************
   PARLib-Method:
   Task:
   Programing:
   07/2004 OK Implementation
   05/2006 WW Fix bugs and add the case of DOF>1
   09/2007 WW Improve the extremely slow performance of the old code
**************************************************************************/
	void CPARDomain::CreateElements()
	{
		//----------------------------------------------------------------------
		if(!m_msh)
			return;
		//----------------------------------------------------------------------
		long i,k;
		int j, nNodes, nNodesHQ;
		long* elem_nodes = NULL;
		MeshLib::IElement* m_ele = NULL;
		MeshLib::INode* m_nod = NULL;
		//*** Buffer for acceleration. 14.09.2007 WW:
		// As long buffer
        std::vector<int> node_connected_doms;
		node_connected_doms.resize((long)m_msh->getNumberOfNodes());
		for(k = 0; k < (long)node_connected_doms.size(); k++)
			node_connected_doms[k] = -1;
		for(k = 0; k < (long)nodes.size(); k++)
		{
			i = nodes[k];
			node_connected_doms[i] = k;
		}
		//***
		//----------------------------------------------------------------------
		for(i = 0; i < (long)elements.size(); i++)
		{
			if(elements[i] > (long)m_msh->getNumberOfElements())
			{
				cout << "Warning: no ELE data" << '\n';
				continue;
			}
			m_ele = m_msh->getElemenet(i);
			nNodes = m_ele->getNumberOfNodes(1); //WW
			nNodesHQ = m_ele->getNumberOfNodes(use_quad?2:1);
			// cout << i << " " << elements[i] << ": ";
			elem_nodes = new long[nNodesHQ]; //WW
			element_nodes_dom.push_back(elem_nodes); //WW
			for(j = 0; j < nNodes; j++)
			{
				elem_nodes[j] = node_connected_doms[m_ele->getNodeID(j)];
			}
			//------------------WW
			if(!use_quad) 
                continue;
			for(j = nNodes; j < nNodesHQ; j++)
			{
				k = m_ele->getNodeID(j);
				if(node_connected_doms[k] > -1)
					elem_nodes[j] = node_connected_doms[k];
				else
				{
					elem_nodes[j] = (long)nodes.size();
					node_connected_doms[k] = elem_nodes[j];
					nodes.push_back(k);
				}
			}
			nnodesHQ_dom = (long) nodes.size();
			//------------------WW
			// cout << endl;
		}
		//
		//----------------------------------------------------------------------
	}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2006 WW Implementation
   09/2007 WW Improve computational efficiency
**************************************************************************/
	void CPARDomain::NodeConnectedNodes()
	{
#if 0
		vector<long> nodes2node;
		// node_connected_doms as buffer to accelerate the computation
		// 14.09.2007 WW
		const size_t n_mesh_nodes (m_msh->nod_vector.size());
		for (size_t i = 0; i < n_mesh_nodes; i++)
			node_connected_doms[i] = -1;

		for (size_t j = 0; j < nodes.size(); j++)
			node_connected_doms[m_msh->nod_vector[nodes[j]]->GetIndex()] = j;

		const size_t n_nodes(nodes.size());
		for (size_t i = 0; i < n_nodes; i++)
		{
			nodes2node.clear();
			std::vector<size_t> const& connected_nodes(
			        m_msh->nod_vector[nodes[i]]->getConnectedNodes());
			const size_t n_connected_nodes(connected_nodes.size());
			for (size_t k = 0; k < n_connected_nodes; k++)
			{
				int j = node_connected_doms[connected_nodes[k]];
				if (j > -1)
					nodes2node.push_back(j);
			}
			const size_t i_buff (nodes2node.size());
			long* nodes_to_node = new long[i_buff];
			for (size_t k = 0; k < i_buff; k++)
				nodes_to_node[k] = nodes2node[k];
			node_conneted_nodes.push_back(nodes_to_node);
			num_nodes2_node.push_back(i_buff);
		}
#endif
	}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   07/2004 OK Implementation
**************************************************************************/
long CPARDomain::GetDOMNode(long global_node)
{
    long i;
    long no_nodes = (long)nodes.size();
    for(i = 0; i < no_nodes; i++)
        if(nodes[i] == global_node)
            return i;
    return -1;
}

/**************************************************************************
   PARLib-Method:
   Task:  Create local equation system
   Programing:
   12/2007 WW Implementation
**************************************************************************/
void CPARDomain::CreateEQS()
{
    if (use_linear) {
        _vec_sparse.push_back(new SparseTable(*this, false));
    }
    if (use_quad) {
        _vec_sparse.push_back(new SparseTable(*this, true));
    }

    _problem2eqs.resize(_problems.size());

    for (size_t i=0; i<_problems.size(); i++) {
        ITransientProblem* p = _problems[i];
        size_t order = 1;
        size_t dof = 1;
        std::pair<size_t, size_t> eqs_type = make_pair(order, dof);
        size_t eqs_id = 0;
        if (_set_eqs.count(eqs_type)==0) {
            Linear_EQS* eqs = new Linear_EQS(*_vec_sparse[order-1], dof);
            _vec_eqs.push_back(eqs);
            eqs_id = _vec_eqs.size()-1;
        } else {
            eqs_id = _set_eqs[eqs_type];
        }
        _problem2eqs[i] = eqs_id;
    }
 }

/**************************************************************************
   PARLib-Method:
   Task:
   Programing:
   12/2007 WW Implementation
**************************************************************************/
void CPARDomain::InitialEQS(size_t problem_id)
{
	Linear_EQS* this_eqs = _vec_eqs[_problem2eqs[problem_id]];
	
	this_eqs->Initialize();
}

//

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   07/2004 OK Implementation
**************************************************************************/
/* //WW
   void CPARDomain::AssembleMatrix(CRFProcess* m_pcs)
   {
   long i;
   //----------------------------------------------------------------------
   SetZeroLinearSolver(eqs);
   //MXDumpGLS("AssembleMatrix1.txt",1,eqs->b,eqs->x);
   //----------------------------------------------------------------------
   long no_elements = (long)elements.size();
   for(i=0;i<no_elements;i++){
    // virtual function PCSAssembleMatrix(i)
   //MakeElementEntryEQS_ASM(elements[i]->global_number,eqs->b,NULL,this);
   //WW    MakeElementEntryEQS_ASM(i,eqs->b,NULL,this,m_pcs);
   //MXDumpGLS("AssembleMatrix1.txt",1,eqs->b,eqs->x);
   }
   }
 */
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   07/2004 OK Implementation
**************************************************************************/
	bool NodeExists(long node,vector<long>node_vector)
	{
		long i;
		long no_nodes = (long)node_vector.size();
		for(i = 0; i < no_nodes; i++)
			if(node == node_vector[i])
				return true;
		return false;
	}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   09/2005 OK MSH
   last modification:
**************************************************************************/
	void CPARDomain::WriteTecplot(string msh_name)
	{
#if 0
		long i;
		string element_type;
		//----------------------------------------------------------------------
		// GSP
		CGSProject* m_gsp = NULL;
		m_gsp = GSPGetMember("msh");
		if(!m_gsp)
			return;
		//--------------------------------------------------------------------
		// file handling
		char name[10];
		sprintf(name,"%i",ID);
		string dom_name = "DOMAIN";
		dom_name += name;
		string dom_file_name = m_gsp->path + dom_name + TEC_FILE_EXTENSION;
		fstream dom_file (dom_file_name.data(),ios::trunc | ios::out);
		dom_file.setf(ios::scientific,ios::floatfield);
		dom_file.precision(12);
		//--------------------------------------------------------------------
		// MSH
		CFEMesh* m_msh = NULL;
		MeshLib::CElem* m_ele = NULL;
		m_msh = FEMGet(msh_name);
		if(!m_msh)
			return;
		//--------------------------------------------------------------------
		if (!dom_file.good())
			return;
		dom_file.seekg(0L,ios::beg);
		//--------------------------------------------------------------------
		for(i = 0; i < (long)elements.size(); i++)
		{
			m_ele = m_msh->ele_vector[elements[i]];
			if(!m_ele)
				continue;
			switch(m_ele->GetElementType())
			{
			case MshElemType::LINE:
				element_type = "ET = QUADRILATERAL";
				break;
			case MshElemType::QUAD:
				element_type = "ET = QUADRILATERAL";
				break;
			case MshElemType::HEXAHEDRON:
				element_type = "ET = BRICK";
				break;
			case MshElemType::TRIANGLE:
				element_type = "ET = TRIANGLE";
				break;
			case MshElemType::TETRAHEDRON:
				element_type = "ET = TETRAHEDRON";
				break;
			case MshElemType::PRISM:
				element_type = "ET = BRICK";
				break;
			default:
				std::cerr << "CPARDomain::WriteTecplot MshElemType not handled" <<
				std::endl;
			}
		}
		//--------------------------------------------------------------------
		dom_file << "VARIABLES = X,Y,Z,DOM" << endl;
		long no_nodes = (long)m_msh->nod_vector.size();
		dom_file << "ZONE T = " << dom_name << ", " \
		         << "N = " << no_nodes << ", " \
		         << "E = " << (long)elements.size() << ", " \
		         << "F = FEPOINT" << ", " << element_type << endl;
		//......................................................................
		for(i = 0; i < no_nodes; i++)
		{
			double const* const coords (m_msh->nod_vector[i]->getData());
			dom_file << coords[0] << " " << coords[1] << " " << coords[2] << " " <<
			ID << endl;
//      m_nod = m_msh->nod_vector[i];
//      dom_file << m_nod->X() << " " << m_nod->Y() << " " << m_nod->Z() << " " << ID << endl;
		}
		//......................................................................
		for(i = 0; i < (long)elements.size(); i++)
		{
			m_ele = m_msh->ele_vector[elements[i]];
			if(!m_ele)
				continue;
			switch(m_ele->GetElementType())
			{
			case MshElemType::LINE:
				dom_file \
				<< m_ele->getNodeIndices()[0] + 1 << " " << m_ele->getNodeIndices()[1] + 1 <<
				" " << m_ele->getNodeIndices()[1] + 1 << " " << m_ele->getNodeIndices()[0] +
				1 << endl;
				element_type = "ET = QUADRILATERAL";
				break;
			case MshElemType::QUAD:
				dom_file \
				<< m_ele->getNodeIndices()[0] + 1 << " " << m_ele->getNodeIndices()[1] + 1 <<
				" " << m_ele->getNodeIndices()[2] + 1 << " " << m_ele->getNodeIndices()[3] +
				1 << endl;
				element_type = "ET = QUADRILATERAL";
				break;
			case MshElemType::HEXAHEDRON:
				dom_file \
				<< m_ele->getNodeIndices()[0] + 1 << " " << m_ele->getNodeIndices()[1] + 1 <<
				" " << m_ele->getNodeIndices()[2] + 1 << " " << m_ele->getNodeIndices()[3] +
				1 << " " \
				<< m_ele->getNodeIndices()[4] + 1 << " " << m_ele->getNodeIndices()[5] + 1 <<
				" " << m_ele->getNodeIndices()[6] + 1 << " " << m_ele->getNodeIndices()[7] +
				1 << endl;
				element_type = "ET = BRICK";
				break;
			case MshElemType::TRIANGLE:
				dom_file \
				<< m_ele->getNodeIndices()[0] + 1 << " " << m_ele->getNodeIndices()[1] + 1 <<
				" " << m_ele->getNodeIndices()[2] + 1 << endl;
				element_type = "ET = TRIANGLE";
				break;
			case MshElemType::TETRAHEDRON:
				dom_file \
				<< m_ele->getNodeIndices()[0] + 1 << " " << m_ele->getNodeIndices()[1] + 1 <<
				" " << m_ele->getNodeIndices()[2] + 1 << " " << m_ele->getNodeIndices()[3] +
				1 << endl;
				element_type = "ET = TETRAHEDRON";
				break;
			case MshElemType::PRISM:
				dom_file \
				<< m_ele->getNodeIndices()[0] + 1 << " " << m_ele->getNodeIndices()[0] + 1 <<
				" " << m_ele->getNodeIndices()[1] + 1 << " " << m_ele->getNodeIndices()[2] +
				1 << " " \
				<< m_ele->getNodeIndices()[3] + 1 << " " << m_ele->getNodeIndices()[3] + 1 <<
				" " << m_ele->getNodeIndices()[4] + 1 << " " << m_ele->getNodeIndices()[5] +
				1 << endl;
				element_type = "ET = BRICK";
				break;
			default:
				std::cerr << "CPARDomain::WriteTecplot MshElemType not handled" <<
				std::endl;
			}
		}
#endif
	}

#if 0
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   10/2005 OK Implementation
   last modification:
**************************************************************************/
	void DOMWriteTecplot(string msh_name)
	{
		CPARDomain* m_dom = NULL;
		for(int i = 0; i < (int)dom_vector.size(); i++)
		{
			m_dom = dom_vector[i];
			m_dom->WriteTecplot(msh_name);
		}
	}
#endif



#if defined(USE_MPI)                              //WW
//------------------------For parallel solvers------------------------------
/*************************************************************************
   GeoSys-Function:
   Task:
   Programming:
   02/2008 WW Implementation
 **************************************************************************/
	void CPARDomain::FillBorderNodeConnectDom(vector<int> allnodes_doms)
	{
		long i, ig;
		int k;
		b_start[0] = num_inner_nodes;
		b_end[0] = num_inner_nodes + num_boundary_nodes;
		nq = 1;
		//
		if(nnodesHQ_dom > nnodes_dom)
		{
			//
			nq = 2;
			b_start[1] = b_end[0] + num_inner_nodesHQ;
			b_end[1] = nnodesHQ_dom;
		}
		//
		//
		for(k = 0; k < nq; k++)
			for(i = b_start[k]; i < b_end[k]; i++)
			{
				ig = nodes[i];
				bnode_connected_dom.push_back(allnodes_doms[ig]);
			}

		/*
		   string test = "rank";
		   static char stro[102];

		   sprintf(stro, "%d",myrank);
		   string test1 = test+(string)stro+"dom.txt";
		   ofstream Dum(test1.c_str(), ios::out); // WW
		   Dum<<b_start[0]<<endl;
		   Dum<<b_end[0]<<endl;
		   Dum<<b_start[1]<<endl;
		   Dum<<b_end[1]<<endl;
		   for(i=0;i<bnode_connected_dom.size();i++)
		   Dum<<bnode_connected_dom[i]<<endl;
		   exit(1);
		 */
	}

/*************************************************************************
   GeoSys-Function:
   Task:
   Programming:
   12/2007 WW Implementation
 **************************************************************************/
	void CPARDomain::ConfigEQS(CNumerics* m_num, const long n, bool quad)
	{
		int i;
		long dim = 0;
		i_start[0] = 0;
		i_end[0] = num_inner_nodes; //Number of interior nodes
		b_start[0] = 0;
		b_end[0] = num_boundary_nodes;
		n_shift[0] = num_inner_nodes;
		long inner_size = num_inner_nodes;
		long border_size = num_boundary_nodes;
		quadratic = quad;
		//
		double cpu_time_local = -MPI_Wtime();

		/*

		   //TEST_MPI
		   string test = "rank";
		   static char stro[102];

		   sprintf(stro, "%d",myrank);
		   string test1 = test+(string)stro+"Assemble.txt";
		   ofstream Dum(test1.c_str(), ios::out); // WW
		   Dum<<"Mysize" <<mysize<<"  "<<quadratic<<endl;

		   Dum.close();
		   MPI_Finalize();
		   exit(1);
		 */

		if(quadratic)
		{
			n_loc = nnodesHQ_dom;
			nq = 2;
			n_bc = t_border_nodes_sizeH;
			i_start[1] = i_end[0] + num_boundary_nodes;
			i_end[1] = i_start[1] + num_inner_nodesHQ; //Number of interior nodes
			//
			b_start[1] = b_end[0];
			b_end[1] = b_start[1] + num_boundary_nodesHQ;
			n_shift[1] =  n_shift[0] + num_inner_nodesHQ;
			//
			dof = eqsH->DOF();
			eqsH->SetDomain(this);
			eqsH->ConfigNumerics(m_num, n);
			inner_size += num_inner_nodesHQ;
			border_size += num_boundary_nodesHQ;
		}
		else
		{
			n_loc = nnodes_dom;
			nq = 1;
			n_bc = t_border_nodes_size;
			dof = eqs->DOF();
			eqs->SetDomain(this);
			eqs->ConfigNumerics(m_num, n);
		}
		//  Concatenate index
		inner_size *= dof;
		border_size *= dof;
		for(i = 0; i < mysize; i++)
		{
			receive_cnt[i] = 1;
			receive_disp[i] = i;
		}
		//
		// receive_cnt_i[]: number of subdomain inner nodes in the concatenated array
		MPI_Allgatherv ( &inner_size, 1, MPI_INT, receive_cnt_i, receive_cnt, receive_disp,
		                 MPI_INT, MPI_COMM_WORLD );
		inner_size = 0;
		for(i = 0; i < mysize; i++)
		{
			receive_disp_i[i] = inner_size;
			inner_size += receive_cnt_i[i];
		}
#if defined(NEW_BREDUCE)
		// receive_cnt_b[]: number of subdomain border nodes in the concatenated array
		MPI_Allgatherv ( &border_size, 1, MPI_INT, receive_cnt_b, receive_cnt, receive_disp,
		                 MPI_INT, MPI_COMM_WORLD );
		border_size = 0;
		for(i = 0; i < mysize; i++)
		{
			receive_disp_b[i] = border_size;
			border_size += receive_cnt_b[i];
		}
		//
		cpu_time_local += MPI_Wtime();
		if(quadratic)
		{
			eqsH->f_buffer[(int)eqsH->f_buffer.size() - 3] = new double[border_size];
			eqsH->cpu_time += cpu_time_local;
		}
		else
		{
			eqs->f_buffer[(int)eqs->f_buffer.size() - 3] = new double[border_size];
			eqs->cpu_time += cpu_time_local;
		}
#endif
		//
		// dim = n_loc*dof;
		// MPI_Allreduce(&dim,  &max_dimen, 1, MPI_INT,  MPI_MAX, MPI_COMM_WORLD);
	}

/*************************************************************************
   GeoSys-Function:
   Task: for parallel solver
   HM monolithic case is to be considered.
   Programming:
   02/2008 WW Implementation
 **************************************************************************/
	double CPARDomain::Dot_Border_Vec(const double* vec_x, const double* vec_y)
	{
		long i, l_buff;
		int ii, k;
		long b_shift[2];
		double val, fac;
		val = 0.;
		//
		for(k = 0; k < nq; k++)
			for(i = b_start[k]; i < b_end[k]; i++)
			{
				fac =  1.0 / (double)bnode_connected_dom[i];
				for(ii = 0; ii < dof; ii++)
				{
					l_buff = i + n_loc * ii + n_shift[k];
					val += fac * vec_x[l_buff] * vec_y[l_buff];
				}
			}
		return val;
	}

/*************************************************************************
   GeoSys-Function:
   Task: Parallel solver

   HM monolithic case is to be considered.
   Programming:
   07/2006 WW Implementation
   12/2007 WW Revise
 **************************************************************************/
	double CPARDomain::Dot_Interior(const double* localr0,  const double* localr1)
	{
		long i;
		int ii, k;
		double val;
		//

		/*
		   if(dof>1)
		   {
		   //TEST

		   string test = "rank";
		   char stro[64];
		   sprintf(stro, "%d",myrank);
		   string test1 = test+(string)stro+"dom.txt";

		   ofstream Dum(test1.c_str(), ios::out);
		   Dum<<" nnodesHQ_dom  "<< nnodesHQ_dom<<endl;

		   Dum<<" nq "<<nq <<endl;

		   for(k=0; k<nq; k++)
		   {
		   Dum<<" i_start[k]  "<<i_start[k] <<endl;

		   Dum<<"  i_end[k] "<<i_end[k] <<endl;

		   for(i=i_start[k];i<i_end[k];i++)
		   {
		   for(ii=0; ii<dof; ii++)
		   {
		   //
		   val += localr0[i+n_loc*ii]*localr0[i+n_loc*ii];

		   Dum<<"[i+n_loc*ii] "<< i+n_loc*ii <<" localr0[i+n_loc*ii] "<< localr0[i+n_loc*ii]<<endl;

		   }
		   }
		   }
		   exit(1);

		   }
		 */

		val = 0.0;
		if(!localr1)
			for(k = 0; k < nq; k++)
				for(i = i_start[k]; i < i_end[k]; i++)
					for(ii = 0; ii < dof; ii++)
						//
						val +=
						        localr0[i + n_loc *
						                ii] * localr0[i + n_loc * ii];

		else
			for(k = 0; k < nq; k++)
				for(i = i_start[k]; i < i_end[k]; i++)
					for(ii = 0; ii < dof; ii++)
						val +=
						        localr0[i + n_loc *
						                ii] * localr1[i + n_loc * ii];

		/*
		   //TEST
		   Dum.close();
		   if(nq>1)
		   exit(1);
		 */

		return val;
	}

/*************************************************************************
   GeoSys-Function:
   Task: Parallel solver
    n: Dimension of the global EQS
   HM monolithic case is to be considered.
   Programming:
   06/2006 WW Implementation
   12/2007 WW Revise
 **************************************************************************/
	void CPARDomain::Global2Local(const double* global_x, double* local_x, const long n )
	{
		long i, ig;
		int ii;
		//
		//
		long n_global = (long)n / dof;
		for(i = 0; i < n_loc * dof; i++)
			local_x[i] = 0.;
		for(i = 0; i < n_loc; i++)
		{
			ig = nodes[i];
			for(ii = 0; ii < dof; ii++)
				local_x[i + n_loc * ii] = global_x[ig + n_global * ii];
		}
	}

/*************************************************************************
   GeoSys-Function:
   Task: Parallel solver
      n: Dimension of the global EQS
   HM monolithic case is to be considered.
   Programming:
   06/2006 WW Implementation
   12/2007 WW Revise
   02/2008 WW Revise
 **************************************************************************/
	void CPARDomain:: Local2Global(const double* local_x, double* global_x, const long n )
	{
		long i, ig, b_index;
		int ii, k;
		double fac = 0.;
		//
		//
		long n_global = (long)n / dof;
		//
		for(i = 0; i < n; i++)
			global_x[i] = 0.;
		//
		for(k = 0; k < nq; k++)
			for(i = i_start[k]; i < i_end[k]; i++)
			{
				ig = nodes[i];
				for(ii = 0; ii < dof; ii++)
					global_x[ig + n_global * ii] = local_x[i + n_loc * ii];
			}
		//
		for(k = 0; k < nq; k++)
			for(i = b_start[k]; i < b_end[k]; i++)
			{
				b_index = i + n_shift[k];
				fac =  1.0 / (double)bnode_connected_dom[i];
				//
				ig = nodes[b_index];
				for(ii = 0; ii < dof; ii++)
					global_x[ig + n_global *
					         ii] += fac * local_x[b_index + n_loc * ii];
			}
		//
	}

/*************************************************************************
   GeoSys-Function:
   Task: Parallel solver

   HM monolithic case is to be considered.
   Programming:
   12/2007 WW
 **************************************************************************/
	void CPARDomain::Global2Border(const double* x, double* local_x, const long n )
	{
		long i, ig;
		int ii, k;
		//
		long nnodes_g = (long)n / dof;
		// BC
		for(i = 0; i < dof * n_bc; i++)
			local_x[i] = 0.0;
		//
		for(i = 0; i < n_bc; i++)
		{
			k = t_border_nodes[i];
			for(ii = 0; ii < dof; ii++)
				local_x[i + n_bc * ii] = x[k + nnodes_g * ii];
		}
	}

/*************************************************************************
   GeoSys-Function:
   Task: Parallel solver

   HM monolithic case is to be considered.
   Programming:
   12/2007 WW
 **************************************************************************/
	void CPARDomain::Border2Global(const double* local_x, double* x, const long n)
	{
		long i, ig;
		int ii, k;
		//
		long nnodes_g = (long)n / dof;
		//
		for(i = 0; i < n_bc; i++)
		{
			k = t_border_nodes[i];
			for(ii = 0; ii < dof; ii++)
				x[k + nnodes_g * ii] = local_x[i + n_bc * ii];
		}
	}

/*************************************************************************
   GeoSys-Function:
   Task: Parallel solver

   HM monolithic case is to be considered.
   Programming:
   07/2006 WW Implementation
   12/2007 WW Revise
 **************************************************************************/
	void CPARDomain::Local2Border(const double* local_x, double* border_x)
	{
		long i, ig;
		int ii, k;
		//
		// BC
		for(i = 0; i < dof * n_bc; i++)
			border_x[i] = 0.0;
		//
		for(k = 0; k < nq; k++)
			for(i = b_start[k]; i < b_end[k]; i++)
			{
				ig = nodes_halo[i];
				for(ii = 0; ii < dof; ii++)
					border_x[ig + n_bc *
					         ii] = local_x[i + n_shift[k] + n_loc * ii];
			}
	}

/*************************************************************************
   GeoSys-Function:
   Task: Parallel solver

   HM monolithic case is to be considered.
   Programming:
   07/2006 WW Implementation
   12/2007 WW Revise
 **************************************************************************/
	void CPARDomain::Border2Local(const double* border_x, double* local_x)
	{
		long i, ig;
		int ii, k;
		//
		//
		for(k = 0; k < nq; k++)
			for(i = b_start[k]; i < b_end[k]; i++)
			{
				ig = nodes_halo[i];
				for(ii = 0; ii < dof; ii++)
					local_x[i + n_shift[k] + n_loc *
					        ii] = border_x[ig + n_bc * ii];
			}
	}

/*\!
 ********************************************************************
   Concatenate the inertanal entries of local subdomain solution
   Programm:
   12/2007 WW
 ********************************************************************/
//#define NEW_BREDUCE2
	void CPARDomain::CatInnerX(double* global_x, const double* local_x, const long n)
	{
		long i, j, ig;
		int ii, k;
		long counter = 0;
		double* x_i, * x_g;
		Linear_EQS* eq = NULL;
		if(quadratic)
			eq = eqsH;
		else
			eq = eqs;
		//
		x_i = eq->f_buffer[0];
		x_g = eq->f_buffer[(long)eq->f_buffer.size() - 1];
		//
		long n_global = (long)n / dof;
		//
#if defined(NEW_BREDUCE2)
		// Not finished
		// Due to the parallel computing of dom topology, not all num_inner_node of
		//   dom_vector[j] is caculated.
		// for(i=0; i<eq->A->Dim();i++)
		//   x_i[i] = 0.;
		//
		for(k = 0; k < nq; k++)
		{
			for(i = i_start[k]; i < i_end[k]; i++)
				for(ii = 0; ii < dof; ii++)
				{
					x_i[counter] = local_x[i + n_loc * ii];
					counter++; //
				}
		}
		// Concatentate
		MPI_Allgatherv (x_i, counter, MPI_DOUBLE, x_g, receive_cnt_i, receive_disp_i,
		                MPI_DOUBLE, MPI_COMM_WORLD );
		//
		// Mapping to the golbal x
		CPARDomain* a_dom;
		for(j = 0; j < mysize; j++)
		{
			counter = receive_disp_i[j];
			// Problem from here
			if(j == myrank)
				a_dom = this;
			else
				;
			a_dom = dom_vector[j];
			a_dom->nq = 1;
			a_dom->i_start[0] = 0;
			a_dom->i_end[0] =  a_dom->num_inner_nodes; //Number of interior nodes
			if(quadratic)
			{
				a_dom->nq = 2;
				a_dom->i_start[1] = a_dom->i_end[0] + a_dom->num_boundary_nodes;
				a_dom->i_end[1] = a_dom->i_start[1] + a_dom->num_inner_nodesHQ;
			}

			for(k = 0; k < a_dom->nq; k++) // This should come from different processors

				for(i = a_dom->i_start[k]; i < a_dom->i_end[k]; i++)
				{
					ig = a_dom->nodes[i];
					for(ii = 0; ii < dof; ii++)
					{
						global_x[ig + n_global * ii] = x_g[counter];
						counter++;
					}
				}
		}

#else // if defined(NEW_BREDUCE2)
		Local2Global(local_x, x_g, n);
		MPI_Allreduce( x_g, global_x, n, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
	}

#if defined(NEW_BREDUCE)
/*\!
 ********************************************************************
   Reduce border entries by concatenating
   Programm:
   12/2007 WW
 ********************************************************************/
	void CPARDomain::ReduceBorderV(double* local_x)
	{
		long i, j, ig;
		int ii, k;
		long counter = 0;
		double* x_b, * x_cat, * x_g;
		Linear_EQS* eq = NULL;
		if(quadratic)
			eq = eqsH;
		else
			eq = eqs;
		//
		x_g = eq->f_buffer[(long)eq->f_buffer.size() - 1];
		x_cat = eq->f_buffer[(int)eq->f_buffer.size() - 2];
		x_b = &local_x[eq->Dim()];
		//
		for(k = 0; k < nq; k++)
		{
			for(i = b_start[k]; i < b_end[k]; i++)
				for(ii = 0; ii < dof; ii++)
				{
					//;
					x_g[counter] = local_x[i + n_shift[k] + n_loc * ii];
					counter++;
				}
		}
		// Concatentate
		MPI_Allgatherv (x_g, counter, MPI_DOUBLE, x_cat, receive_cnt_b, receive_disp_b,
		                MPI_DOUBLE, MPI_COMM_WORLD );
		for(i = 0; i < dof * n_bc; i++)
			x_b[i] = 0.0;
		//
		CPARDomain* a_dom;
		for(j = 0; j < mysize; j++)
		{
			a_dom = dom_vector[j];
			counter = receive_disp_b[j];
			for(k = 0; k < a_dom->nq; k++)
				for(i = a_dom->b_start[k]; i < a_dom->b_end[k]; i++)
				{
					ig = a_dom->nodes_halo[i];
					for(ii = 0; ii < dof; ii++)
					{
						x_b[ig + n_bc * ii] += x_cat[counter];
						counter++;
					}
				}
		}
	}
#endif                                            //if defined(NEW_BREDUCE)
/********************************************************************
   As the title
   Programm:
   12/2007 WW
********************************************************************/
	void CPARDomain::PrintEQS_CPUtime(ostream &os)
	{
		if(eqs)
			os << "CPU time elapsed in linear solver for linear elements: "
			   << eqs->GetCPUtime() << endl;
		if(eqsH)
			os << "CPU time elapsed in linear solver for quadratic elements: "
			   << eqsH->GetCPUtime() << endl;
	}
#endif                                            //// if defined(USE_MPI)
}
}
