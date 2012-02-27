
#include "par_ddc_group.h"

#include <iostream>

#define MAX_ZEILE 512
#define DDC_FILE_EXTENSION ".ddc"

using namespace std;

namespace OGS5
{

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   07/2004 OK Implementation
   07/2004 OK Version 1 untill 3.9.17OK6
   07/2004 OK Version 2 from   3.9.17OK7
   10/2005 OK cout
**************************************************************************/
void CPARDomainGroup::DOMRead(string file_base_name)
{
	//----------------------------------------------------------------------
	std::cout << "DOMRead: ";
	//----------------------------------------------------------------------
	CPARDomain* m_dom = NULL;
	char line[MAX_ZEILE];
	string sub_line;
	string line_string;
	string ddc_file_name;
	ios::pos_type position;
	//----------------------------------------------------------------------
	// File handling
	ddc_file_name = file_base_name + DDC_FILE_EXTENSION;
	ifstream ddc_file (ddc_file_name.data(),ios::in);
	if(!ddc_file.good())
	{
		cout << "no DDC file" << endl;
		return;
	}
	ddc_file.seekg(0L,ios::beg);
	//----------------------------------------------------------------------
	// Keyword loop
	while (!ddc_file.eof())
	{
		ddc_file.getline(line,MAX_ZEILE);
		line_string = line;
		//----------------------------------------------------------------------
		// keyword found
		if(line_string.find("#DOMAIN") != string::npos)
		{
			m_dom = new CPARDomain();
			position = m_dom->Read(&ddc_file);
			dom_vector.push_back(m_dom);
			ddc_file.seekg(position,ios::beg);
		}                         // keyword found
	}                                     // eof
	//----------------------------------------------------------------------
	cout << dom_vector.size() << " domains" << endl;
	//----------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2007 WW Implementation
**************************************************************************/
void CPARDomainGroup::CountDoms2Nodes(bool quad)
{
	size_t i,j,k;
	long n_index = 0;
	CPARDomain* m_dom = NULL;
	MeshLib::INode* anode = NULL;
	MeshLib::IElement* elem = NULL;
	MeshLib::IMesh* a_msh = m_msh;

	//Average of nodal Neumann BCs contributed by nodes from different domains
	size_t nsize = a_msh->getNumberOfNodes();
	node_connected_doms.resize(nsize);

	for(j = 0; j < nsize; j++)
		node_connected_doms[j] = 0;
    std::vector<bool> nod_mark(nsize, false);
    std::vector<bool> ele_mark(a_msh->getNumberOfElements(), false);
	for(i = 0; i < dom_vector.size(); i++) {
		m_dom = dom_vector[i];
		for(j = 0; j < m_dom->elements.size(); j++) {
			elem = a_msh->getElemenet(m_dom->elements[j]);
			if(ele_mark[elem->getID()]) {
                for(k = 0; k < elem->getNumberOfNodes(quad); k++) {
                    n_index = elem->getNodeID(k);
                    if(!nod_mark[n_index]) {
                        node_connected_doms[n_index] += 1;
                        nod_mark[n_index] = true;
                    }
                }
            }
		}
	}
	//
	for(j = 0; j < a_msh->getNumberOfElements(); j++)
	{
		elem = a_msh->getElemenet(j);
		if(ele_mark[j])
			for(k = 0; k < elem->getNumberOfNodes(quad); k++) 
			{
                nod_mark[elem->getNodeID(k)] = true;
			}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   07/2004 OK Implementation
   05/2006 WW Fix bugs and add the case of DOF>1
   07/2006 WW Find nodes of all neighbors of each node
   09/2007 WW Check the nodes whether it is belong to the deactivated elements
   10/2010 TF changed access to process type
**************************************************************************/

void CPARDomainGroup::setup()
{
	const size_t no_domains = dom_vector.size();
	if(no_domains == 0) return;

	CPARDomain* m_dom = NULL;
	size_t i;

	//----------------------------------------------------------------------
	// Create domain nodes
	cout << "->Create DOM" << endl;
    for(i = 0; i < no_domains; i++)
    {
        m_dom = dom_vector[i];
        m_dom->setMesh(this->m_msh, use_linear, use_quad);
//        m_dom->setProblems(_problems);
    }

	//----------------------------------------------------------------------
	// Create domain nodes
	cout << "  Create domain nodes" << endl;
	for(i = 0; i < no_domains; i++)
	{
		m_dom = dom_vector[i];
		cout << "    Domain:" << m_dom->ID << endl;
		m_dom->CreateNodes();
	}
	//----------------------------------------------------------------------
	// Create domain elements
	cout << "  Create domain elements" << endl;
	for(i = 0; i < no_domains; i++)
	{
		m_dom = dom_vector[i];
		cout << "    Domain:" << m_dom->ID << endl;
		m_dom->CreateElements();
	}
	// For find nodes connected to node WW
	long j;
	long nsize = m_msh->getNumberOfNodes();
	//node_connected_doms.resize(nsize);
	for(j = 0; j < nsize; j++)
		node_connected_doms[j] = 0;
	for(i = 0; i < dom_vector.size(); i++)
	{
		m_dom = dom_vector[i];
		for(j = 0; j < (long)m_dom->nodes.size(); j++)
			node_connected_doms[m_dom->nodes[j]] += 1;
	}
	//
	// Find nodes of all neighbors of each node. // WW
	// Local topology. WW
	cout << "  Find nodes on borders" << endl;
	FindNodesOnInterface(use_quad);
	cout << "  Find the connected nodes for each node" << endl;
#ifndef USE_MPI                                //WW
	for(i = 0; i < no_domains; i++)
	{
#else
	i = myrank;                           //WW
#endif
		m_dom = dom_vector[i];
		m_dom->NodeConnectedNodes();
#ifndef USE_MPI                             //WW
	}
#endif
	//----------------------------------------------------------------------
	// Create domain EQS
	cout << "  Create domain EQS" << endl;
#ifdef USE_MPI
    i = myrank;
#else
	for(i = 0; i < no_domains; i++)
	{
#endif
		m_dom = dom_vector[i];
		cout << "    Domain:" << m_dom->ID << endl;
		m_dom->CreateEQS();
		//
#ifndef USE_MPI
	}
#endif
		//----------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   07/2006 WW Implementation
   last modification:
**************************************************************************/
	void CPARDomainGroup::FindNodesOnInterface(bool quadr)
	{
		// int k;
		long i, nnodes_gl, g_index, nnodes_l;
		size_t j;
		long l_buff = 0, l_buff1 = 0;

		long* elem_nodes = NULL;
		//
        MeshLib::IElement* m_ele = 0;
		//
		CPARDomain* m_dom = NULL;
		vector<long> boundary_nodes;
		vector<long> boundary_nodes_HQ;
		vector<long> inner_nodes_HQ;
		vector<long> inner_nodes;
		vector<long> long_buffer;
		vector<long> bc_buffer;
		vector<long> dom_bc_buffer;
		vector<long> dom_bc_bufferHQ;
		//
		nnodes_gl = m_msh->getNumberOfNodes(); //TODO total nodes
		nnodes_l = m_msh->getNumberOfNodes();	//
		bc_buffer.resize(nnodes_gl);
		//
#if defined(USE_MPI)                           //13.12.2007
		long overlapped_entry_size = 0;
#endif
		for(i = 0; i < nnodes_gl; i++)
		{
			if(node_connected_doms[i] > 1.0)
			{
				// Mapp BC entry-->global node array
				bc_buffer[i] = (long)long_buffer.size();
				long_buffer.push_back(i);
#ifdef USE_MPI
				if(m_msh->nod_vector[i]->GetIndex() < m_msh->GetNodesNumber(false))
					overlapped_entry_size = (long)long_buffer.size();
#endif
			}
			else
				bc_buffer[i] = -1;
		}
		//
#if defined(USE_MPI)                           // 13.12.2007
		// Total border nodes
		m_dom = dom_vector[myrank];
		m_dom->t_border_nodes_size = overlapped_entry_size;
		m_dom->t_border_nodes_sizeH = (long)long_buffer.size();
		m_dom->t_border_nodes = new long[m_dom->t_border_nodes_sizeH];
		for(i = 0; i < m_dom->t_border_nodes_sizeH; i++)
			m_dom->t_border_nodes [i] = long_buffer[i];
#endif

		// Sort
#ifndef USE_MPI
		for(int k = 0; k < (int)dom_vector.size(); k++)
		{
			int myrank = k;
#endif
		m_dom = dom_vector[myrank];
		//
		boundary_nodes.clear();
		inner_nodes.clear();
		boundary_nodes_HQ.clear();
		inner_nodes_HQ.clear();
		long_buffer.clear();
		long_buffer.resize((long)m_dom->nodes.size());
		dom_bc_buffer.clear();
		dom_bc_bufferHQ.clear();

		for(i = 0; i < (long)m_dom->nodes.size(); i++)
		{
			g_index = m_dom->nodes[i];
			if(node_connected_doms[ g_index] > 1.0)
			{
				if(g_index >= nnodes_l)
				{
					boundary_nodes_HQ.push_back(i);
					dom_bc_bufferHQ.push_back(bc_buffer[g_index]);
					long_buffer[i] = -(long)boundary_nodes_HQ.size() -
					                 nnodes_gl;
				}
				else
				{
					boundary_nodes.push_back(i);
					dom_bc_buffer.push_back(bc_buffer[g_index]);
					long_buffer[i] = -(long)boundary_nodes.size();
				}
			}
			else
			{
				if(g_index >= nnodes_l)
				{
					long_buffer[i] = (long)inner_nodes_HQ.size() + nnodes_gl;
					inner_nodes_HQ.push_back(i);
				}
				else
				{
					long_buffer[i] = (long)inner_nodes.size();
					inner_nodes.push_back(i);
				}
			}
		}
		//
		m_dom->num_inner_nodes = (long)inner_nodes.size();
		m_dom->num_inner_nodesHQ = (long)inner_nodes_HQ.size();
		m_dom->num_boundary_nodes = (long)boundary_nodes.size();
		m_dom->num_boundary_nodesHQ = (long)boundary_nodes_HQ.size();
		// Sort for high order nodes

		m_dom->nodes_inner.clear();
		m_dom->nodes_halo.clear();
		for(i = 0; i < m_dom->num_inner_nodes; i++)
			m_dom->nodes_inner.push_back(m_dom->nodes[inner_nodes[i]]);
		for(i = 0; i < (long)inner_nodes_HQ.size(); i++)
			m_dom->nodes_inner.push_back(m_dom->nodes[inner_nodes_HQ[i]]);
		//
		for(i = 0; i < m_dom->num_boundary_nodes; i++)
			m_dom->nodes_halo.push_back(m_dom->nodes[boundary_nodes[i]]);
		for(i = 0; i < (long)boundary_nodes_HQ.size(); i++)
			m_dom->nodes_halo.push_back(m_dom->nodes[boundary_nodes_HQ[i]]);
		//
		m_dom->nodes.clear();
		m_dom->nodes.resize(m_dom->nnodesHQ_dom);
		// First interior nodes, then interface nodes
		j = 0;
		for(i = 0; i < m_dom->num_inner_nodes; i++)
			m_dom->nodes[i] = m_dom->nodes_inner[i];
		//       m_dom->nodes_inner[i] = i;
		j += m_dom->num_inner_nodes;
		for(i = 0; i < m_dom->num_boundary_nodes; i++)
			m_dom->nodes[i + j] = m_dom->nodes_halo[i];
		//      m_dom->nodes_halo[i] = i+j;
		j += m_dom->num_boundary_nodes;
		for(i = 0; i < (long)inner_nodes_HQ.size(); i++)
			m_dom->nodes[i + j] = m_dom->nodes_inner[i + m_dom->num_inner_nodes];
		//     m_dom->nodes_inner[i+m_dom->num_inner_nodes] = i+j;
		j += (long)inner_nodes_HQ.size();
		for(i = 0; i < (long)boundary_nodes_HQ.size(); i++)
			m_dom->nodes[i + j] = m_dom->nodes_halo[i + m_dom->num_boundary_nodes];
		//      m_dom->nodes_halo[i+m_dom->num_boundary_nodes] = i+j;

		for(i = 0; i < (long)m_dom->nodes.size(); i++)
		{
			l_buff = long_buffer[i];
			if(l_buff < 0) //interface nodes
			{
				if(-l_buff - nnodes_gl > 0) //HQ nodes
					l_buff1 = m_dom->nnodes_dom + (long)inner_nodes_HQ.size() -
					          l_buff - nnodes_gl - 1;
				else
					l_buff1 = (long)inner_nodes.size() - l_buff - 1;
				//             l_buff1 = m_dom->num_inner_nodesHQ-l_buff-1;
			}
			else
			{
				if(l_buff - nnodes_gl >= 0) //HQ nodes
					l_buff1 = m_dom->nnodes_dom + l_buff - nnodes_gl;
				else
					l_buff1 = l_buff;
			}
			long_buffer[i] = l_buff1;
		}
		//
#ifdef USE_MPI                              //WW
		m_dom->FillBorderNodeConnectDom(node_connected_doms);
#endif
		//
		m_dom->nodes_inner.clear();
		m_dom->nodes_halo.clear();
		// Mapping the local index to global BC array, overlapped_entry.
		//
		for(i = 0; i < m_dom->num_boundary_nodes; i++)
			m_dom->nodes_halo.push_back(dom_bc_buffer[i]);
		for(i = 0; i < (long)boundary_nodes_HQ.size(); i++)
			m_dom->nodes_halo.push_back(dom_bc_bufferHQ[i]);
		//----------------------------------------------------------------------
		for(i = 0; i < (long)m_dom->elements.size(); i++)
		{
			m_ele = m_msh->getElemenet(m_dom->elements[i]);
			elem_nodes = m_dom->element_nodes_dom[i];
			for(j = 0; j < m_ele->getNumberOfNodes(quadr); j++)
			{
				l_buff = elem_nodes[j];
				elem_nodes[j] = long_buffer[l_buff];
			}
		}

#ifndef USE_MPI
	}
#endif
	}

/**************************************************************************
   DDCLib-Function
   07/2007 OK Encapsulation
   10/2010 TF changed access to process type
 ***************************************************************************/
	void DDCCreate()
	{
#if 0
		//----------------------------------------------------------------------
		// DDC
		if(dom_vector.size() > 0)
		{
			//WW ----- Domain decomposition ------------------
			int i;
			int no_processes = (int)pcs_vector.size();
			CRFProcess* m_pcs = NULL;
			bool DOF_gt_one = false;
			//----------------------------------------------------------------------
			for(i = 0; i < no_processes; i++)
			{
				m_pcs = pcs_vector[i];
				//if(m_pcs->pcs_type_name.find("DEFORMATION")!=string::npos) { // TF 10/2010
				if(m_pcs->getProcessType () == FiniteElement::DEFORMATION ||
				   m_pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW)
				{
					DOF_gt_one = true;
					break;
				}
			}
			if(!DOF_gt_one)
				m_pcs = pcs_vector[0];
			// -----------------------
			DOMCreate();
			//
			for(i = 0; i < no_processes; i++)
			{
				m_pcs = pcs_vector[i];
				// Config boundary conditions for domain decomposition
				m_pcs->SetBoundaryConditionSubDomain(); //WW
			}
			//
			node_connected_doms.clear();
		}
		// PA PCSProcessDependencies();
#endif
	}


} //end