#pragma once

#include "SingleFlowDCRBendersSolver.h"
#include "ChangeSolver.h"
#include <queue>

namespace SMSpp_di_unipi_it
{

    class SFDCRBChange : public Change
    {
    public:
        // valutare, dato che so se è fissata o meno
        enum SFDCRBChangeType
        {
            eFixX = 0,                // fix an arc to not be used in the solution
            eUnfixX = 1,              // restore an arc to the solution
            eLastSFDCRBChangeType = 2 // first allowed new change type for derived classes
        };
        // constructor
        SFDCRBChange() = default;

        SFDCRBChange(Block::Index arc, SFDCRBChangeType type) : Change(), f_Arc(arc), f_Type(type) {}

        // destructor
        virtual ~SFDCRBChange() = default;

        // getter

        Block::Index get_arc(void) const
        {
            return f_Arc;
        }

        SFDCRBChangeType get_type(void) const
        {
            return f_Type;
        }

        void serialize(netCDF::NcGroup &group) const override
        {
            Change::serialize(group);
            group.addVar("arc", netCDF::NcType::nc_UINT64);
            group.addVar("type", netCDF::NcType::nc_INT);
            group.getVar("arc").putVar(&f_Arc);
            group.getVar("type").putVar(&f_Type);
        }

        void deserialize(const netCDF::NcGroup &group) override
        {
            group.getVar("arc").getVar(&f_Arc);
            group.getVar("type").getVar(&f_Type);
        }

        Change *apply(Block *block, bool doUndo = false, ModParam issueMod = eNoBlck, ModParam issueAMod = eNoBlck) override
        {
            auto B = dynamic_cast<SingleFlowDCRBlock *>(block);
            if (!B)
                throw(std::invalid_argument("SFDCRBChange::apply: block is not a SingleFlowDCRBlock"));
            if (f_Arc >= B->get_NArcs())
                throw(std::invalid_argument("SFDCRBChange::apply: arc index is out of bounds"));

            if (f_Type == eFixX && B->is_closed(f_Arc))
            {
                throw std::invalid_argument("SFDCRBChange::apply: arc is already fixed");
            }
            else if (f_Type == eUnfixX && !B->is_closed(f_Arc))
            {
                throw std::invalid_argument("SFDCRBChange::apply: arc is already unfixed");
            }
            else if (f_Type == eFixX)
            {
                B->close_arc(f_Arc, issueMod, issueAMod);
            }
            else if (f_Type == eUnfixX)
            {
                B->open_arc(f_Arc, issueMod, issueAMod);
            }
            else
            {
                throw std::runtime_error("Unknown SFDCRBChangeType");
            }

            if (doUndo)
            {
                return new SFDCRBChange(f_Arc, f_Type == eFixX ? eUnfixX : eFixX);
            }
            else
            {
                return nullptr;
            }
        }

    protected:
        Block::Index f_Arc;      // the arc to be fixed/unfixed
        SFDCRBChangeType f_Type; // the type of change (fix/unfix)

    private:
        SMSpp_insert_in_factory_h;
    };

    // Gives access to the protected state of DCRLagrangianSolver, nothing more
    class SFDCRBProbe : public DCRLagrangianSolver
    {
    public:
        /// after Solve() at a fixed r_min: the two paths that define lambda*,
        /// as arc indices of the *full* arc set (0-based, source-to-sink
        /// order); what is found is written on log, if any
        bool two_paths(std::vector<int> &Pplus, std::vector<int> &Pminus,
                       std::ostream *log = nullptr)
        {
            Solve();

            if (log)
                *log << "[probe] status=" << lagstat << " lambda=" << lambda
                     << " pLambda=" << pLambda << " mLambda=" << mLambda
                     << " rmin=" << rmin << std::endl;

            if ((lagstat != DCRLagrangianSolver::OK) || (lambda <= 0))
                return (false);

            const double ls = lambda; // path_at() overwrites lambda
            const double d = 1e-6 * ls;

            Pminus = path_at(ls - d); // just below: P- (beta >= 0)
            Pplus = path_at(ls + d);  // just above: P+ (beta <  0)

            if (Pminus == Pplus)
            { // fallback: the two cut points
                Pminus = path_at(pLambda);
                Pplus = path_at(mLambda);
            }

            if (log)
                *log << "[probe] |P-|=" << Pminus.size()
                     << " |P+|=" << Pplus.size()
                     << " equal=" << (Pminus == Pplus) << std::endl;

            return ((!Pminus.empty()) && (!Pplus.empty()) &&
                    (Pminus != Pplus));
        }
        /// true if, at the loaded r_min, lambda = 0 is optimal and its path feasible
        bool zero_opt(void)
        {
            if (isFeasible()) // returns 0 when the relaxation is OK
                return (false);
            return ((getLambda() == 0) && (getHeurVal() < Inf<double>()));
        }

    private:
        std::vector<int> path_at(double lam)
        {
            lambda = lam;
            setSPTcosts(); // same costs as when the cut was generated
            spt.updCosts(Linksp);
            spt.Solve();
            if (spt.getStatus() != SPT::OK)
                return (std::vector<int>());

            auto p = spt.getPath();
            std::vector<int> out(p.size());
            for (std::size_t k = 0; k < p.size(); ++k)
                out[k] = RedGraPos[p[k]]; // reduced index -> arc index
            return (out);
        }
    };

    /// the RelaxationSolver of the nodes of a BranchAndXSolver on a
    /// SingleFlowDCRBlock; when a log stream is set [see Solver::set_log()],
    /// each compute(), branch() and certificate of the node writes one line
    /// on it, and nothing is written otherwise
    class SFDCRBRelaxationSolver : public SingleFlowDCRBendersSolver,
                                   public RelaxationSolver
    {

    public:
        // constructor
        SFDCRBRelaxationSolver() : SingleFlowDCRBendersSolver() {}

        // destructor
        virtual ~SFDCRBRelaxationSolver() = default;

        /*--------------------------------------------------------------------------*/
        void set_Block(Block *block) override
        {
            SingleFlowDCRBendersSolver::set_Block(block);
            f_closed.clear();
            f_infeasible = false;
            if (auto B = dynamic_cast<SingleFlowDCRBlock *>(block))
                f_closed.assign(B->get_NArcs(), 0);
        }

        /*--------------------------------------------------------------------------*/

        // COMPUTE: reload BenBound with the arcs this node has closed, then solve

        int compute(bool changedvars = true) override
        {
            if (!f_Block)
                return (kBlockLocked);

            const bool owned = f_Block->is_owned_by(f_id);
            if ((!owned) && (!f_Block->read_lock()))
                return (kBlockLocked);

            // Block changes first: the base would reload them *without* our closures
            process_outstanding_Modification();
            load_node(); // always: BenBound::Solve() cannot be re-run
                         // on the same data (solvedflag), it needs a
                         // fresh LoadProblem()
            if (f_log)
                *f_log << "[compute] rho=" << f_flow.rate
                       << " burst=" << f_flow.burst
                       << " deadline=" << f_flow.deadline << " infeasible="
                       << f_infeasible << std::endl;
            if (!owned)
                f_Block->read_unlock();

            if (f_infeasible) // no s-t path left: nothing to solve
                return (kOK);
            const int rc = SingleFlowDCRBendersSolver::compute(changedvars);

            f_cert = -SMSpp_di_unipi_it::Inf<OFValue>();
            if ((!f_infeasible) &&
                (get_ub() - SingleFlowDCRBendersSolver::get_lb() >
                 1e-6 * std::max(1.0, std::abs(get_ub()))))
                certify_left();

            return (rc);
        }

        /*--------------------------------------------------------------------------*/
        // APPLY: closures live in the solver, the Block is not touched

        Change *apply(Change *chg, bool doUndo = false) override
        {
            if (auto ac = dynamic_cast<SFDCRBChange *>(chg))
            {
                const auto a = ac->get_arc();
                if (a >= f_closed.size())
                    throw(std::invalid_argument("SFDCRBRelaxationSolver::apply: bad arc"));

                const bool was_closed = f_closed[a]; // state BEFORE acting
                f_closed[a] = (ac->get_type() == SFDCRBChange::eFixX);

                return (doUndo ? new SFDCRBChange(a, was_closed ? SFDCRBChange::eFixX
                                                                : SFDCRBChange::eUnfixX)
                               : nullptr);
            }

            return (chg->apply(f_Block, doUndo)); // anything else: the Block
        }

        /*--------------------------------------------------------------------------*/
        // BRANCH: first arc where the feasible and the infeasible path split

        std::vector<Change *> branch(void) override
        {
            if (f_log)
                *f_log << "[branch] stat=" << BenBound::getStat()
                       << " r*=" << BenBound::getr_min()
                       << " LB=" << get_lb() << " UB=" << get_ub()
                       << " heur=" << BenBound::getHeurVal()
                       << " benIt=" << BenBound::getNumIterationBender()
                       << " narcs=" << f_narcs << std::endl;
            if (f_infeasible || (BenBound::getStat() != BenBound::OK))
                throw(std::logic_error("branch(): the node has no relaxation"));

            SFDCRBProbe pr;
            pr.LoadProblem(f_nnodes, f_narcs, f_flow, f_links.data(),
                           f_nodes.data(), f_mtu, BenBound::getr_min());

            std::vector<int> Pp, Pm;
            if (!pr.two_paths(Pp, Pm, f_log))
                throw(std::logic_error("branch(): no pair of distinct paths at r*"));

            std::size_t k = 0; // common prefix, arc by arc
            while ((k < Pp.size()) && (k < Pm.size()) && (Pp[k] == Pm[k]))
                ++k;
            if ((k >= Pp.size()) || (k >= Pm.size()))
                throw(std::logic_error("branch(): one path is a prefix of the other"));

            // both arcs leave the same node u, and a simple path uses at most one
            // of them: every routing avoids a+ or avoids a-, so the children cover
            // the node. Closing a- first keeps the feasible path P+ in the first
            // child (a heuristic order: it favours the incumbent, test it)
            return (std::vector<Change *>{new SFDCRBChange(Pm[k], SFDCRBChange::eFixX),
                                          new SFDCRBChange(Pp[k], SFDCRBChange::eFixX)});
        }

        /*--------------------------------------------------------------------------*/

        // getter
        OFValue get_lb(void) override
        {
            if (f_infeasible || (BenBound::getStat() != BenBound::OK))
                return (SMSpp_di_unipi_it::Inf<OFValue>());
            return (std::max(SingleFlowDCRBendersSolver::get_lb(), f_cert));
        }
        OFValue get_true_lb(void) override { return (get_lb()); }
        OFValue get_true_ub(void) override { return (get_ub()); }
        bool has_true_var_solution(void) override { return (has_var_solution()); }

        void get_true_var_solution(Configuration *solc = nullptr) override
        {
            get_var_solution(solc);
        }

        Solution *get_true_solution(Configuration *solc = nullptr) override
        {
            return (get_Solution(solc));
        }

    private:
        void certify_left(void)
        {
            f_cert = -SMSpp_di_unipi_it::Inf<OFValue>();
            if (f_infeasible)
                return;

            const double rlow = std::max(f_flow.rate,
                                         f_flow.burst / f_flow.deadline);
            SFDCRBProbe pr;
            pr.LoadProblem(f_nnodes, f_narcs, f_flow, f_links.data(),
                           f_nodes.data(), f_mtu, rlow);
            if (pr.zero_opt())
                f_cert = pr.getObjVal(); // optimum of the node
            if (f_log)
                *f_log << "[cert] rlow=" << rlow << " f_cert=" << f_cert
                       << std::endl;
        }

        void load_node(void)
        {
            auto DCRB = static_cast<SingleFlowDCRBlock *>(f_Block);
            f_nnodes = DCRB->get_NNodes();
            f_narcs = DCRB->get_NArcs();
            f_mtu = DCRB->get_MTU();
            if (f_closed.size() != f_narcs)
                f_closed.assign(f_narcs, 0);

            int source = -1, sink = -1;
            for (Index i = 0; i < f_nnodes; ++i)
            {
                const auto Bi = DCRB->get_B(i);
                if (Bi < 0)
                    source = i;
                else if (Bi > 0)
                    sink = i;
            }
            if ((source < 0) || (sink < 0))
                throw(std::logic_error("load_node: no source or no sink"));

            f_flow = {};
            f_flow.sourcenode = source;
            f_flow.sinknode = sink;
            f_flow.burst = DCRB->get_FlowBurst();
            f_flow.rate = DCRB->get_rho();
            f_flow.deadline = DCRB->get_FlowDeadline();

            // no feasible r_min is below this, so an arc with less capacity than
            // half of it can never be in the reduced graph
            const double rlow = std::max(f_flow.rate,
                                         f_flow.burst / f_flow.deadline);
            if (!(rlow > 0))
                throw(std::logic_error("load_node: needs rho > 0 or burst > 0"));

            auto &LD = DCRB->get_LinkDelays();
            auto &ND = DCRB->get_NodeDelays();

            f_links.assign(f_narcs, {});
            std::vector<std::vector<Index>> adj(f_nnodes);
            for (Index i = 0; i < f_narcs; ++i)
            {
                const bool closed = f_closed[i] || DCRB->is_closed(i) ||
                                    DCRB->is_deleted(i);
                auto &l = f_links[i];
                l.startnode = DCRB->get_SN(i) - 1;
                l.endnode = DCRB->get_EN(i) - 1;
                l.speed = DCRB->get_U(i);
                l.capacity = closed ? 0.5 * rlow : DCRB->get_U(i);
                l.delay = LD.empty() ? 0 : LD[i];
                l.cost = DCRB->is_deleted(i) ? 0 : DCRB->get_C(i);
                if (!closed)
                    adj[l.startnode].push_back(l.endnode);
            }

            f_nodes.assign(f_nnodes, {});
            for (Index i = 0; i < f_nnodes; ++i)
                f_nodes[i].delay = ND.empty() ? 0 : ND[i];

            // is the sink still reachable on the open arcs?
            std::vector<char> seen(f_nnodes, 0);
            std::queue<Index> q;
            q.push(source);
            seen[source] = 1;
            while (!q.empty())
            {
                const auto u = q.front();
                q.pop();
                for (auto v : adj[u])
                    if (!seen[v])
                    {
                        seen[v] = 1;
                        q.push(v);
                    }
            }
            f_infeasible = !seen[sink];

            BenBound::LoadProblem(f_nnodes, f_narcs, f_flow, f_links, f_nodes,
                                  f_mtu);
        }

        std::vector<char> f_closed; ///< arcs closed by apply()
        bool f_infeasible = false;  ///< no s-t path on the open arcs
        DCR::DCRFlow f_flow;
        std::vector<DCR::DCRLink> f_links;
        std::vector<DCR::DCRNode> f_nodes;
        Index f_nnodes = 0, f_narcs = 0;
        double f_mtu = 0;
        OFValue f_cert = -SMSpp_di_unipi_it::Inf<OFValue>(); ///< lower bound at r_min = rlow

        SMSpp_insert_in_factory_h;
    };
}