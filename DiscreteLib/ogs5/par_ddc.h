
#ifndef par_dd_INC
#define par_dd_INC


// C++ STL
#include <fstream>
#include <list>
#include <string>
#include <vector>
#include <map>

#include "MeshLib/Core/IMesh.h"


namespace OGS5
{

class SparseTable;
class Linear_EQS;

/**
 * \brief Sub-domain
 */
class CPARDomain
{
public:
    CPARDomain(size_t id, MeshLib::IMixedOrderMesh &msh, bool linear, bool quad);
    ~CPARDomain(void);

    /// get domain id
    size_t getID() const {return _dom_id;};
    /// get mesh
    MeshLib::IMixedOrderMesh* getMesh() const {return _msh;};

    /// get the number of this domain nodes
    long getNumberOfDomainNodes(bool quad) const { return quad ? _nnodesHQ_dom : _nnodes_dom; }
    /// get the total number of this domain nodes
    size_t getTotalNumberOfDomainNodes() const {return _list_dom_global_nodes.size();};
    /// resize the number domain nodes
    void resizeDomainNodes(size_t n) {_list_dom_global_nodes.resize(n);};
    /// get local node id
    long getLocalNodeID(size_t global_id) const;
    /// get global node id
    long getGlobalNodeID(size_t local_id) const {return _list_dom_global_nodes[local_id];};
    /// set global node id
    void setGlobalNodeID(size_t local_id, long global_id) {_list_dom_global_nodes[local_id] = global_id;};

    /// get the number of inner nodes
    long getNumberOfInnerNodes(bool quad) const {return quad ? _num_inner_nodesHQ : _num_inner_nodes;};
    /// set the number of inner nodes
    void setNumberOfInnerNodes(bool quad, long n) { if (quad) _num_inner_nodesHQ = n; else _num_inner_nodes = n;};
    /// reset a list of inner nodes
    void resetInnerNodes() { _list_inner_nodes_global.clear(); }
    /// add inner node
    void addInnerNode(long global_id) {_list_inner_nodes_global.push_back(global_id);};
    /// get inner node
    long getInnerNode(size_t i) const {return _list_inner_nodes_global[i];};

    /// get the number of boundary nodes
    long getNumberOfBoundaryNodes(bool quad) const {return quad ? num_boundary_nodesHQ : num_boundary_nodes;};
    /// set the number of boundary nodes
    void setNumberOfBoundaryNodes(bool quad, long n) { if (quad) num_boundary_nodesHQ = n; else num_boundary_nodes = n;};
    /// reset a list of boundary nodes
    void resetBoundaryNodes() { _list_boundary_nodes_global.clear(); }
    /// add boundary node
    void addBoundaryNode(long global_id) {_list_boundary_nodes_global.push_back(global_id);};
    /// get boundary node
    long getBoundaryNode(size_t i) const {return _list_boundary_nodes_global[i];};

    /// get the number of domain elements
    size_t getNumberOfElements() const {return _list_dom_elements.size(); };
    /// get element id
    long getElementId(size_t i) const {return _list_dom_elements[i];};

    /// get the number of connected nodes 
    long getNumberOfNodesConnectedToNode(size_t i) const
    {
        return _node2conneted_nodes[i].size();
    }
    /// set the number of connected nodes
    void setNumberOfNodesConnectedToNode(size_t i, long v)
    {
        _node2conneted_nodes[i].resize(v);
    }
    long get_node_conneted_nodes(size_t i, size_t j) const
    {
        return _node2conneted_nodes[i][j];
    }
    void set_node_conneted_nodes(size_t i, size_t j, long v)
    {
        _node2conneted_nodes[i][j] = v;
    }

    long* get_element_nodes_dom(size_t i) const {return element_nodes_dom[i];};

    void NodeConnectedNodes();


    void CreateEQS();
	void InitialEQS(size_t problem_id);
	void CalcElementMatrices();

    Linear_EQS* getEQS(bool quad)
    {
        if (_vec_eqs.size()==0) return 0;
        if (!quad) return _vec_eqs[0];
        if (_use_linear) return _vec_eqs[1];
        else return _vec_eqs[0];
    }

#if defined(USE_MPI)                           //WW
	// long MaxDim() const {return max_dimen;}   //WW
	void ReleaseMemory();
	//WW
	void FillBorderNodeConnectDom(std::vector<int> allnodes_doms);
	long BSize() const                    //WW
	{
		return n_bc;
	}
	//WW
	void ConfigEQS(CNumerics* m_num, const long n, bool quad = false);
	//WW
	double Dot_Interior(const double* localr0,  const double* localr1 = NULL);
	//WW
	void Global2Local(const double* global_x, double* local_x, const long n );
	//WW
	void Local2Global(const double* local_x, double* global_x, const long n );
	//WW
	void Global2Border(const double* x, double* local_x, const long n);
	//WW
	void Border2Global(const double* local_x, double* x, const long n);
	//WW
	void Local2Border(const double* local_x, double* border_x);
	//WW
	void Border2Local(const double* border_x, double* local_x);
	//
	double Dot_Border_Vec(const double* vec_x, const double* vec_y);
	//
	//WW
	void CatInnerX(double* global_x, const double* local_x, const long n);
	//WW
	void PrintEQS_CPUtime(std::ostream &os = std::cout);

#if defined(NEW_BREDUCE)
	void ReduceBorderV(double* local_x);
#endif
	//
	//
#endif

private:
    void CreateNodes();
    void CreateElements();

    std::map<std::pair<size_t, size_t>, size_t> _set_eqs;
    std::vector<Linear_EQS*> _vec_eqs;
    std::vector<SparseTable*> _vec_sparse;

    std::vector<long*> element_nodes_dom; // Local DOM element nodes. WW
    long _nnodes_dom;
    long _nnodesHQ_dom;
    //#ifdef USE_MPI //WW
    long _num_inner_nodes;
    long _num_inner_nodesHQ;
    long num_boundary_nodes;
    long num_boundary_nodesHQ;
    //#endif
#if defined(USE_MPI)                           // 13.12.2007 WW
    // Store global indices of all border nodes to border_nodes of the whole mesh
    // 0-->border_nodes_size, nodes for linear interpolation
    // border_nodes_size-->border_nodes_sizeH, nodes for quadratic interpolation
    std::vector<int> bnode_connected_dom; // Connected doms of border nodes
    long* t_border_nodes;
    long t_border_nodes_size;
    long t_border_nodes_sizeH;
    // For local EQS
    // For index mapping from local to global
    int dof, nq;
    long n_loc, n_bc;                     //, max_dimen;
    long i_start[2], i_end[2];
    long b_start[2], b_end[2], n_shift[2];
    // Data for concatenate internal entries of domain vectors
    int* receive_cnt_i;
    int* receive_disp_i;
#if defined(NEW_BREDUCE)
    // Data for concatenate border entries of domain vectors
    int* receive_cnt_b;
    int* receive_disp_b;
#endif
    //
    int* receive_cnt;
    int* receive_disp;
    //friend class Math_Group::Linear_EQS;
    //
#endif
    //
    // Equation
    int _dom_id;
    std::vector<long> _list_dom_elements;
    std::vector<long> _list_inner_nodes_global;
    std::vector<long> _list_boundary_nodes_global;
    std::vector<long> _list_dom_global_nodes;
    MeshLib::IMixedOrderMesh* _msh;
    bool _quadratic;
    std::vector<std::vector<long>> _node2conneted_nodes;
    bool _use_linear;
    bool _use_quad;
};

}

#endif
