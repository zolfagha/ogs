
#include "par_ddc_group.h"

#include <iostream>

#define MAX_ZEILE 512
#define DDC_FILE_EXTENSION ".ddc"

using namespace std;

namespace OGS5
{


void CPARDomainGroup::countDoms2Nodes(bool quad)
{
	//Average of nodal Neumann BCs contributed by nodes from different domains
	const size_t nsize = _msh->getNumberOfNodes();
	_node_connected_doms.resize(nsize, 0);

    std::vector<bool> nod_mark(nsize, false);
    std::vector<bool> ele_mark(_msh->getNumberOfElements(), false);
	for (size_t i=0; i<_dom_vector.size(); i++) {
		CPARDomain* m_dom = _dom_vector[i];
		for (size_t j=0; j<m_dom->getNumberOfElements(); j++) {
			MeshLib::IElement* elem = _msh->getElemenet(m_dom->getElementId(j));
            if (!ele_mark[elem->getID()]) continue;
            for(size_t k=0; k<elem->getNumberOfNodes(quad); k++) {
                size_t n_index = elem->getNodeID(k);
                if (!nod_mark[n_index]) {
                    _node_connected_doms[n_index] += 1;
                    nod_mark[n_index] = true;
                }
            }
		}
	}
	
	for (size_t j=0; j<_msh->getNumberOfElements(); j++) {
		MeshLib::IElement* elem = _msh->getElemenet(j);
		if(ele_mark[j]) {
			for (size_t k=0; k<elem->getNumberOfNodes(quad); k++) {
                nod_mark[elem->getNodeID(k)] = true;
			}
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
	const size_t no_domains = _dom_vector.size();
	if(no_domains == 0) return;

	CPARDomain* m_dom = NULL;
	size_t i;

	// For find nodes connected to node WW
	long j;
	long nsize = _msh->getNumberOfNodes();
	//node_connected_doms.resize(nsize);
	for(j = 0; j < nsize; j++)
		_node_connected_doms[j] = 0;
	for(i = 0; i < _dom_vector.size(); i++)
	{
		m_dom = _dom_vector[i];
		for(j = 0; j < (long)m_dom->getTotalNumberOfDomainNodes(); j++)
			_node_connected_doms[m_dom->getGlobalNodeID(j)] += 1;
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
		m_dom = _dom_vector[i];
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
		m_dom = _dom_vector[i];
		cout << "    Domain:" << m_dom->getID() << endl;
		m_dom->CreateEQS();
		//
#ifndef USE_MPI
	}
#endif
		//----------------------------------------------------------------------
}

void CPARDomainGroup::findNodesOnInterface(bool quadr)
{
	const size_t nnodes_gl = _msh->getNumberOfTotalNodes(); 
	const size_t nnodes_l = _msh->getNumberOfNodes();	//
	
#if defined(USE_MPI)
	long overlapped_entry_size = 0;
#endif
    // Map BC entry-->global node array
    vector<long> list_bndNodeId;
    vector<long> map_glNode2bndNodeId(nnodes_gl);
	for (size_t i=0; i<nnodes_gl; i++) {
		if (_node_connected_doms[i] > 1) {
			map_glNode2bndNodeId[i] = (long)list_bndNodeId.size();
			list_bndNodeId.push_back(i);
#ifdef USE_MPI
			if (_msh->nod_vector[i]->GetIndex() < _msh->GetNodesNumber(false))
				overlapped_entry_size = (long)list_bndNodeId.size();
#endif
		} else {
			map_glNode2bndNodeId[i] = -1;
        }
	}
	
#if defined(USE_MPI)
	// Total border nodes
	m_dom = _dom_vector[myrank];
	m_dom->t_border_nodes_size = overlapped_entry_size;
	m_dom->t_border_nodes_sizeH = (long)list_bndNodeId.size();
	m_dom->t_border_nodes = new long[m_dom->t_border_nodes_sizeH];
	for(long i = 0; i < m_dom->t_border_nodes_sizeH; i++)
		m_dom->t_border_nodes [i] = list_bndNodeId[i];
#endif

	// Sort
#ifndef USE_MPI
	for (size_t k=0; k<_dom_vector.size(); k++) {
		int myrank = (int) k;
#endif
        CPARDomain* m_dom = _dom_vector[myrank];
        setupDomain(m_dom, map_glNode2bndNodeId, quadr);
#ifndef USE_MPI
    }
#endif
}

void CPARDomainGroup::setupDomain( CPARDomain* m_dom, const vector<long>& map_glNode2bndNodeId, bool quadr ) 
{
    //
    const size_t nnodes_gl = _msh->getNumberOfTotalNodes(); 
    const size_t nnodes_l = _msh->getNumberOfNodes();
    const size_t n_dom_nodes = m_dom->getTotalNumberOfDomainNodes();

    // find inner nodes, boundary nodes
    vector<size_t> list_local_boundary_nodes;
    vector<size_t> list_local_boundary_nodes_HQ;
    vector<long> list_local_inner_nodes_HQ;
    vector<long> list_local_inner_nodes;
    vector<long> map_localBndNode2globalBndNode;
    vector<long> map_localBndNode2globalBndNode_HQ;
    vector<long> map_localNode2sortedNode(n_dom_nodes);

    for (size_t i=0; i<n_dom_nodes; i++) {
        const long g_index = m_dom->getGlobalNodeID(i);
        if (_node_connected_doms[g_index] > 1) {
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
    m_dom->setNumberOfInnerNodes(false, (long)list_local_inner_nodes.size());
    m_dom->setNumberOfInnerNodes(true, (long)list_local_inner_nodes_HQ.size());
    m_dom->resetInnerNodes();
    for(size_t i=0; i<list_local_inner_nodes.size(); i++)
        m_dom->addInnerNode(m_dom->getGlobalNodeID(list_local_inner_nodes[i]));
    for(size_t i=0; i<list_local_inner_nodes_HQ.size(); i++)
        m_dom->addInnerNode(m_dom->getGlobalNodeID(list_local_inner_nodes_HQ[i]));
    // set boundary nodes
    m_dom->setNumberOfBoundaryNodes(false, (long)list_local_boundary_nodes.size());
    m_dom->setNumberOfBoundaryNodes(true, (long)list_local_boundary_nodes_HQ.size());
    m_dom->resetBoundaryNodes();
    for(size_t i = 0; i < list_local_boundary_nodes.size(); i++)
        m_dom->addBoundaryNode(m_dom->getGlobalNodeID(list_local_boundary_nodes[i]));
    for(size_t i = 0; i < list_local_boundary_nodes_HQ.size(); i++)
        m_dom->addBoundaryNode(m_dom->getGlobalNodeID(list_local_boundary_nodes_HQ[i]));

    // set all nodes in the domain
    m_dom->resizeDomainNodes(m_dom->getNumberOfDomainNodes(true));
    // First interior nodes, then interface nodes
    long j = 0;
    for (long i = 0; i < m_dom->getNumberOfInnerNodes(false); i++)
        m_dom->setGlobalNodeID(i, m_dom->getInnerNode(i));
    j += m_dom->getNumberOfInnerNodes(false);
    for(long i = 0; i < m_dom->getNumberOfBoundaryNodes(false); i++)
        m_dom->setGlobalNodeID(i+j, m_dom->getBoundaryNode(i));
    j += m_dom->getNumberOfBoundaryNodes(false);
    for(long i = 0; i < (long)list_local_inner_nodes_HQ.size(); i++)
        m_dom->setGlobalNodeID(i + j, m_dom->getInnerNode(i + m_dom->getNumberOfInnerNodes(false)));
    j += (long)list_local_inner_nodes_HQ.size();
    for(long i = 0; i < (long)list_local_boundary_nodes_HQ.size(); i++)
        m_dom->setGlobalNodeID(i + j, m_dom->getBoundaryNode(i + m_dom->getNumberOfBoundaryNodes(false)));


    // create a list of sorted node id
    for (size_t i=0; i<m_dom->getTotalNumberOfDomainNodes(); i++) {
        const long signed_local_bnd_node_id = map_localNode2sortedNode[i];
        long corrected_local_node_id = 0;
        if (signed_local_bnd_node_id < 0) { 
            //interface nodes
            if (- signed_local_bnd_node_id - nnodes_gl > 0) //HQ nodes
                corrected_local_node_id = m_dom->getNumberOfDomainNodes(false) + (long)list_local_inner_nodes_HQ.size() - signed_local_bnd_node_id - nnodes_gl - 1;
            else
                corrected_local_node_id = (long)list_local_inner_nodes.size() - signed_local_bnd_node_id - 1;
        } else {
            // internal
            if (signed_local_bnd_node_id - nnodes_gl >= 0) //HQ nodes
                corrected_local_node_id = m_dom->getNumberOfDomainNodes(false) + (signed_local_bnd_node_id - nnodes_gl);
            else
                corrected_local_node_id = signed_local_bnd_node_id;
        }
        map_localNode2sortedNode[i] = corrected_local_node_id;
    }

    //
#ifdef USE_MPI
    m_dom->FillBorderNodeConnectDom(_node_connected_doms);
#endif

    //
    m_dom->resetInnerNodes();
    m_dom->resetBoundaryNodes();
    // Mapping the local index to global BC array, overlapped_entry.
    for(long i=0; i<m_dom->getNumberOfBoundaryNodes(false); i++)
        m_dom->addBoundaryNode(map_localBndNode2globalBndNode[i]);
    for(size_t i=0; i<list_local_boundary_nodes_HQ.size(); i++)
        m_dom->addBoundaryNode(map_localBndNode2globalBndNode_HQ[i]);
    //----------------------------------------------------------------------
    for (size_t i = 0; i<m_dom->getNumberOfElements(); i++) {
        MeshLib::IElement* m_ele = _msh->getElemenet(m_dom->getElementId(i));
        long* elem_nodes = m_dom->get_element_nodes_dom(i);
        for (size_t j = 0; j < m_ele->getNumberOfNodes(quadr); j++) {
            long unsorted_id = elem_nodes[j];
            elem_nodes[j] = map_localNode2sortedNode[unsorted_id];
        }
    }
}

} //end


