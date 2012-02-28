
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
	findNodesOnInterface(use_quad);
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

void CPARDomainGroup::findNodesOnInterface(bool quadr)
{
	const size_t nnodes_gl = m_msh->getNumberOfTotalNodes(); 
	const size_t nnodes_l = m_msh->getNumberOfNodes();	//
	
#if defined(USE_MPI)
	long overlapped_entry_size = 0;
#endif
    // Map BC entry-->global node array
    vector<long> list_bndNodeId;
    vector<long> map_glNode2bndNodeId(nnodes_gl);
	for (size_t i=0; i<nnodes_gl; i++) {
		if (node_connected_doms[i] > 1) {
			map_glNode2bndNodeId[i] = (long)list_bndNodeId.size();
			list_bndNodeId.push_back(i);
#ifdef USE_MPI
			if (m_msh->nod_vector[i]->GetIndex() < m_msh->GetNodesNumber(false))
				overlapped_entry_size = (long)list_bndNodeId.size();
#endif
		} else {
			map_glNode2bndNodeId[i] = -1;
        }
	}
	
#if defined(USE_MPI)
	// Total border nodes
	m_dom = dom_vector[myrank];
	m_dom->t_border_nodes_size = overlapped_entry_size;
	m_dom->t_border_nodes_sizeH = (long)list_bndNodeId.size();
	m_dom->t_border_nodes = new long[m_dom->t_border_nodes_sizeH];
	for(long i = 0; i < m_dom->t_border_nodes_sizeH; i++)
		m_dom->t_border_nodes [i] = list_bndNodeId[i];
#endif

	// Sort
#ifndef USE_MPI
	for (size_t k=0; k<dom_vector.size(); k++) {
		int myrank = (int) k;
#endif
        CPARDomain* m_dom = dom_vector[myrank];
        setupDomain(m_dom, map_glNode2bndNodeId, quadr);
#ifndef USE_MPI
    }
#endif
}

void CPARDomainGroup::setupDomain( CPARDomain* m_dom, const vector<long>& map_glNode2bndNodeId, bool quadr ) 
{
    //
    const size_t nnodes_gl = m_msh->getNumberOfTotalNodes(); 
    const size_t nnodes_l = m_msh->getNumberOfNodes();

    // find inner nodes, boundary nodes
    vector<size_t> list_local_boundary_nodes;
    vector<size_t> list_local_boundary_nodes_HQ;
    vector<long> list_local_inner_nodes_HQ;
    vector<long> list_local_inner_nodes;
    vector<long> map_localBndNode2globalBndNode;
    vector<long> map_localBndNode2globalBndNode_HQ;
    vector<long> map_localNode2sortedNode(m_dom->nodes.size());

    for (size_t i=0; i<m_dom->nodes.size(); i++) {
        const long g_index = m_dom->nodes[i];
        if (node_connected_doms[g_index] > 1) {
            // boundary
            if (g_index >= (long)nnodes_l) {
                list_local_boundary_nodes_HQ.push_back(i);
                map_localBndNode2globalBndNode_HQ.push_back(map_glNode2bndNodeId[g_index]);
                map_localNode2sortedNode[i] = -(long)list_local_boundary_nodes_HQ.size() - nnodes_gl;
            } else {
                list_local_boundary_nodes.push_back(i);
                map_localBndNode2globalBndNode.push_back(map_glNode2bndNodeId[g_index]);
                map_localNode2sortedNode[i] = -(long)list_local_boundary_nodes.size();
            }
        } else {
            // internal
            if (g_index >= (long)nnodes_l) {
                map_localNode2sortedNode[i] = (long)list_local_inner_nodes_HQ.size() + nnodes_gl;
                list_local_inner_nodes_HQ.push_back(i);
            } else {
                map_localNode2sortedNode[i] = (long)list_local_inner_nodes.size();
                list_local_inner_nodes.push_back(i);
            }
        }
    }
    // set inner nodes
    m_dom->num_inner_nodes = (long)list_local_inner_nodes.size();
    m_dom->num_inner_nodesHQ = (long)list_local_inner_nodes_HQ.size();
    m_dom->nodes_inner.clear();
    for(long i=0; i<m_dom->num_inner_nodes; i++)
        m_dom->nodes_inner.push_back(m_dom->nodes[list_local_inner_nodes[i]]);
    for(long i=0; i<(long)list_local_inner_nodes_HQ.size(); i++)
        m_dom->nodes_inner.push_back(m_dom->nodes[list_local_inner_nodes_HQ[i]]);
    // set boundary nodes
    m_dom->num_boundary_nodes = (long)list_local_boundary_nodes.size();
    m_dom->num_boundary_nodesHQ = (long)list_local_boundary_nodes_HQ.size();
    m_dom->nodes_halo.clear();
    for(long i = 0; i < m_dom->num_boundary_nodes; i++)
        m_dom->nodes_halo.push_back(m_dom->nodes[list_local_boundary_nodes[i]]);
    for(long i = 0; i < (long)list_local_boundary_nodes_HQ.size(); i++)
        m_dom->nodes_halo.push_back(m_dom->nodes[list_local_boundary_nodes_HQ[i]]);

    // set all nodes in the domain
    m_dom->nodes.clear();
    m_dom->nodes.resize(m_dom->nnodesHQ_dom);
    // First interior nodes, then interface nodes
    long j = 0;
    for (long i = 0; i < m_dom->num_inner_nodes; i++)
        m_dom->nodes[i] = m_dom->nodes_inner[i];
    j += m_dom->num_inner_nodes;
    for(long i = 0; i < m_dom->num_boundary_nodes; i++)
        m_dom->nodes[i + j] = m_dom->nodes_halo[i];
    j += m_dom->num_boundary_nodes;
    for(long i = 0; i < (long)list_local_inner_nodes_HQ.size(); i++)
        m_dom->nodes[i + j] = m_dom->nodes_inner[i + m_dom->num_inner_nodes];
    j += (long)list_local_inner_nodes_HQ.size();
    for(long i = 0; i < (long)list_local_boundary_nodes_HQ.size(); i++)
        m_dom->nodes[i + j] = m_dom->nodes_halo[i + m_dom->num_boundary_nodes];


    // create a list of sorted node id
    for (size_t i=0; i<m_dom->nodes.size(); i++) {
        const long signed_local_bnd_node_id = map_localNode2sortedNode[i];
        long corrected_local_node_id = 0;
        if (signed_local_bnd_node_id < 0) { 
            //interface nodes
            if (- signed_local_bnd_node_id - nnodes_gl > 0) //HQ nodes
                corrected_local_node_id = m_dom->nnodes_dom + (long)list_local_inner_nodes_HQ.size() - signed_local_bnd_node_id - nnodes_gl - 1;
            else
                corrected_local_node_id = (long)list_local_inner_nodes.size() - signed_local_bnd_node_id - 1;
        } else {
            // internal
            if (signed_local_bnd_node_id - nnodes_gl >= 0) //HQ nodes
                corrected_local_node_id = m_dom->nnodes_dom + (signed_local_bnd_node_id - nnodes_gl);
            else
                corrected_local_node_id = signed_local_bnd_node_id;
        }
        map_localNode2sortedNode[i] = corrected_local_node_id;
    }

    //
#ifdef USE_MPI                              //WW
    m_dom->FillBorderNodeConnectDom(node_connected_doms);
#endif
    //
    m_dom->nodes_inner.clear(); //not needed any more?
    m_dom->nodes_halo.clear();
    // Mapping the local index to global BC array, overlapped_entry.
    for(long i=0; i<m_dom->num_boundary_nodes; i++)
        m_dom->nodes_halo.push_back(map_localBndNode2globalBndNode[i]);
    for(size_t i=0; i<list_local_boundary_nodes_HQ.size(); i++)
        m_dom->nodes_halo.push_back(map_localBndNode2globalBndNode_HQ[i]);
    //----------------------------------------------------------------------
    for (size_t i = 0; i<m_dom->elements.size(); i++) {
        MeshLib::IElement* m_ele = m_msh->getElemenet(m_dom->elements[i]);
        long* elem_nodes = m_dom->element_nodes_dom[i];
        for (size_t j = 0; j < m_ele->getNumberOfNodes(quadr); j++) {
            long unsorted_id = elem_nodes[j];
            elem_nodes[j] = map_localNode2sortedNode[unsorted_id];
        }
    }
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


