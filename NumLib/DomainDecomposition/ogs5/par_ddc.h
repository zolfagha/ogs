/**************************************************************************
   PARLib - Object: PAR
   Task: class implementation
   Programing:
   07/2004 OK Implementation
   last modified:
**************************************************************************/
#ifndef par_dd_INC
#define par_dd_INC


// C++ STL
#include <fstream>
#include <list>
#include <string>
#include <vector>
#include <map>

#include "MeshLib/Core/IMesh.h"

#include "NumLib/TimeStepping/ITransientSystem.h"


namespace NumLib
{
namespace OGS5
{

class SparseTable;
class Linear_EQS;

class CPARDomain
{
public:
    CPARDomain(void);
    ~CPARDomain(void);

    std::ios::pos_type Read(std::ifstream*);
    void WriteTecplot(std::string);       //OK

    int getID() const {return ID;};
    void setID(int id) {ID = id;};
    MeshLib::IMixedOrderMesh* getMesh() const {return m_msh;};
    void setMesh(MeshLib::IMixedOrderMesh* msh, bool linear, bool quad) {m_msh = msh; use_linear=linear; use_quad=quad;};
    void setProblems(std::vector<ITransientSystem*> &p) { _problems = p; }


    void CreateNodes();
	void CreateElements();
	void NodeConnectedNodes();            //WW

    void CreateEQS();                     //WW
	void InitialEQS(size_t problem_id);   //WW
	void CalcElementMatrices();

    void assembly(size_t problem_id)
    {
        Linear_EQS *eqs = _vec_eqs[_problem2eqs[problem_id]];
        for (size_t i=0; i<elements.size(); i++) {
            MeshLib::IElement* e = m_msh->getElemenet(i);
            // local assembly
            // add to eqs
        }
    }

    long GetDOMNode(long);
	//
	long GetDomainNodes() const           //WW
	{
		if(quadratic) return nnodesHQ_dom;
		else return nnodes_dom;
	}
	long GetDomainNodes(bool quad) const  //WW
	{
		if(quad) return nnodesHQ_dom;
		else return nnodes_dom;
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
    friend class CPARDomainGroup;
    friend class SparseTable;

    std::vector<ITransientSystem*> _problems;
    std::vector<size_t> _problem2eqs;
    std::map<std::pair<size_t, size_t>, size_t> _set_eqs;
    std::vector<Linear_EQS*> _vec_eqs;
    std::vector<SparseTable*> _vec_sparse;

    std::vector<long*> element_nodes_dom; // Local DOM element nodes. WW
    long nnodes_dom;
    long nnodesHQ_dom;
    //#ifdef USE_MPI //WW
    long num_inner_nodes;
    long num_inner_nodesHQ;
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
    int ID;
    std::vector<long> elements;
    std::vector<long> nodes_inner;
    std::vector<long> nodes_halo;
    std::vector<long> nodes;
    MeshLib::IMixedOrderMesh* m_msh;
    bool selected;                        //OK
    bool quadratic;                       //WW
    std::vector<long*> node_conneted_nodes; //WW
    std::vector<int> num_nodes2_node;     //WW
    int m_color[3];                       //OK
    bool use_linear;
    bool use_quad;
};

}
}

#endif
