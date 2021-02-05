#ifndef NODE_BDD_DUMPER_HPP
#define NODE_BDD_DUMPER_HPP

#include <unordered_set>
#include "NodeBdd.hpp"
#include "util/MyHashTable.hpp"
#include "util/MyList.hpp"
#include "util/MyVector.hpp"

/**
 * DD dumper.
 * A node table is printed in Graphviz (dot) format.
 */
template <typename S, typename T = NodeBdd<double>>
class DdDumper {
    using Spec = S;
    static int const AR = Spec::ARITY;
    static int const headerSize = 1;

    /* SpecNode
     * ┌────────┬────────┬────────┬─────
     * │ nodeId │state[0]│state[1]│ ...
     * └────────┴────────┴────────┴─────
     */
    struct SpecNode {
        NodeId nodeId;
    };

    static NodeId& nodeId(SpecNode* p) { return p->nodeId; }

    static NodeId nodeId(SpecNode const* p) { return p->nodeId; }

    static void* state(SpecNode* p) { return p + headerSize; }

    static void const* state(SpecNode const* p) { return p + headerSize; }

    template <typename SPEC>
    struct Hasher {
        SPEC const& spec;
        int const   level;

        Hasher(SPEC const& _spec, int _level) : spec(_spec), level(_level) {}

        size_t operator()(SpecNode const* p) const {
            return spec.hash_code(state(p), level);
        }

        size_t operator()(SpecNode const* p, SpecNode const* q) const {
            return spec.equal_to(state(p), state(q), level);
        }
    };

    using UniqTable = std::unordered_set<SpecNode*, Hasher<Spec>, Hasher<Spec>>;

    static int getSpecNodeSize(int n) {
        if (n < 0) {
            throw std::runtime_error("storage size is not initialized!!!");
        }
        return headerSize + (n + sizeof(SpecNode) - 1) / sizeof(SpecNode);
    }

    Spec      spec;
    int const specNodeSize;
    char*     oneState;
    NodeId    oneId;

    std::vector<MyList<SpecNode>> spec_nodes_table;
    std::vector<UniqTable>        uniqTable;
    std::vector<Hasher<Spec>>     hasher;

   public:
    explicit DdDumper(Spec const& s)
        : spec(s),
          specNodeSize(getSpecNodeSize(spec.datasize())),
          oneState(nullptr),
          oneId(1) {}

    ~DdDumper() {
        if (oneState) {
            spec.destruct(oneState);
            delete[] oneState;
        }
    }

    DdDumper<S, T>(DdDumper<S, T>&&) = default;
    DdDumper<S, T>(const DdDumper<S, T>&) = default;
    DdDumper<S, T>& operator=(const DdDumper<S, T>&) = default;
    DdDumper<S, T>& operator=(DdDumper<S, T>&&) = default;

    /**
     * Dumps the node table in Graphviz (dot) format.
     * @param os the output stream.
     * @param title title label.
     */
    void dump(std::ostream& os, const std::string& title) {
        if (oneState) {
            spec.destruct(oneState);
        } else {
            oneState = new char[spec.datasize()];
        }
        int n = spec.get_root(oneState);

        os << "digraph \"" << title << "\" {\n";

        if (n == 0) {
            if (!title.empty()) {
                os << "  labelloc=\"t\";\n";
                os << "  label=\"" << title << "\";\n";
            }
        } else if (n < 0) {
            os << R"(  "^" [shape=none,label=")" << title << "\"];\n";
            os << R"(  "^" -> ")" << oneId << "\" [style=dashed"
               << "];\n";
            os << "  \"" << oneId << "\" ";
            os << "[shape=square,label=\"⊤\"];\n";
        } else {
            NodeId root(n, 0);

            for (int i = n; i >= 1; --i) {
                os << "  " << i << " [shape=none,label=\"";
                spec.printLevel(os, i);
                os << "\"];\n";
            }
            for (int i = n - 1; i >= 1; --i) {
                os << "  " << (i + 1) << " -> " << i << " [style=invis];\n";
            }

            os << R"(  "^" [shape=none,label=")" << title << "\"];\n";
            auto dummy = R"(  "^" -> ")";
            os << dummy << root << "\" [style=dashed"
               << "];\n";

            spec_nodes_table.clear();
            spec_nodes_table.resize(n + 1);
            SpecNode* p = spec_nodes_table[n].alloc_front(specNodeSize);
            spec.destruct(oneState);
            spec.get_copy(state(p), oneState);
            nodeId(p) = root;

            uniqTable.clear();
            uniqTable.reserve(n + 1);
            hasher.clear();
            hasher.reserve(n + 1);
            for (int i = 0; i <= n; ++i) {
                hasher.push_back(Hasher<Spec>(spec, i));
                uniqTable.push_back(UniqTable(0, hasher.back(), hasher.back()));
            }

            for (int i = n; i >= 1; --i) {
                dumpStep(os, i);
            }

            for (size_t j = 2; j < oneId.code(); ++j) {
                os << "  \"" << NodeId(j) << "\" ";
                os << "[style=invis];\n";
            }
            os << "  \"" << oneId << "\" ";
            os << "[shape=square,label=\"⊤\"];\n";
        }

        os << "}\n";
        os.flush();
    }

   private:
    void dumpStep(std::ostream& os, int i) {
        MyList<SpecNode>& spec_nodes = spec_nodes_table[i];
        size_t const      m = spec_nodes.size();
        std::vector<char> tmp(spec.datasize());
        void* const       tmpState = tmp.data();
        std::vector<T>    nodeList(m);

        for (size_t j = m - 1; j + 1 > 0; --j, spec_nodes.pop_front()) {
            NodeId f(i, j);
            assert(!spec_nodes.empty());
            SpecNode* p = spec_nodes.front();

            os << "  \"" << f << "\" [label=\"";
            spec.print_state(os, state(p), i);
            os << "\"];\n";

            for (int b = 0; b < AR; ++b) {
                NodeId& child = nodeList[j].branch[b];

                if (nodeId(p) == 0) {
                    child = 0;
                    continue;
                }

                spec.get_copy(tmpState, state(p));
                int ii = spec.get_child(tmpState, i, b);

                if (ii == 0) {
                    child = 0;
                } else if (ii < 0) {
                    if (oneId == 1) {  // the first 1-terminal candidate
                        oneId = 2;
                        spec.destruct(oneState);
                        spec.get_copy(oneState, tmpState);
                        child = oneId;
                    } else {
                        switch (spec.merge_states(oneState, tmpState)) {
                            case 1:
                                oneId = oneId.code() + 1;
                                spec.destruct(oneState);
                                spec.get_copy(oneState, tmpState);
                                child = oneId;
                                break;
                            case 2:
                                child = 0;
                                break;
                            default:
                                child = oneId;
                                break;
                        }
                    }
                } else {
                    SpecNode* pp =
                        spec_nodes_table[ii].alloc_front(specNodeSize);
                    size_t jj = spec_nodes_table[ii].size() - 1;
                    spec.get_copy(state(pp), tmpState);

                    // SpecNode*& pp0 = uniqTable[ii].add(pp);
                    auto aux = uniqTable[ii].insert(pp);
                    if (aux.second) {
                        nodeId(pp) = child = NodeId(ii, jj);
                    } else {
                        auto pp0 = *(aux.first);
                        switch (spec.merge_states(state(pp0), state(pp))) {
                            case 1:
                                nodeId(pp0) = 0;
                                nodeId(pp) = child = NodeId(ii, jj);
                                pp0 = pp;
                                break;
                            case 2:
                                child = 0;
                                spec.destruct(state(pp));
                                spec_nodes_table[ii].pop_front();
                                break;
                            default:
                                child = nodeId(pp0);
                                spec.destruct(state(pp));
                                spec_nodes_table[ii].pop_front();
                                break;
                        }
                    }
                }

                spec.destruct(tmpState);
            }

            spec.destruct(state(p));
        }

        for (size_t j = 0; j < m; ++j) {
            for (int b = 0; b < AR; ++b) {
                NodeId f(i, j);
                NodeId child = nodeList[j].branch[b];
                if (child == 0) {
                    continue;
                }
                if (child == 1) {
                    child = oneId;
                }

                os << "  \"" << f << "\" -> \"" << child << "\"";

                os << " [style=";
                if (b == 0) {
                    os << "dashed";
                } else {
                    os << "solid";
                    if (AR > 2) {
                        os << ",color="
                           << ((b == 1)   ? "blue"
                               : (b == 2) ? "red"
                                          : "green");
                    }
                }
                os << "];\n";
            }
        }

        os << "  {rank=same; " << i;
        for (size_t j = 0; j < m; ++j) {
            os << "; \"" << NodeId(i, j) << "\"";
        }
        os << "}\n";

        uniqTable[i - 1].clear();
        // spec.destructLevel(i);
    }
};

#endif  // NODE_BDD_DUMPER_HPP
