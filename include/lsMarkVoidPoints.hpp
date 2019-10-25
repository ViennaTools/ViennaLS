#ifndef LS_MARK_VOID_POINTS_HPP
#define LS_MARK_VOID_POINTS_HPP

#include <lsPreCompileMacros.hpp>

#include <lsGraph.hpp>

/// This class is used to mark points of the level set
/// which are enclosed in a void.
template <class T, int D> class lsMarkVoidPoints {
  const lsDomain<T, D> *domain = nullptr;

  apply() {
    boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> Graph;

    unsigned int num_components = 0;
    // unsigned int total_number_of_runs=0;

    // allocate memory for component list
    //        std::vector<int> comp_lst[l.number_of_segments()][D+1];
    //        std::vector<int> *comp_lst = new std::vector<int>
    //        [l.number_of_segments()][D+1];


    // std::vector<std::vector<int>> comp_lst(domain->getNumberOfSegments());
    // for(auto it = comp_lst.begin(); it != comp_lst.end(); ++it){
    //   it->resize();
    // }

    std::vector<int> **comp_lst =
        new std::vector<int> *[l.number_of_segments()];
    for (unsigned int i = 0; i < l.number_of_segments(); ++i) {
      comp_lst[i] = new std::vector<int>[D + 1];
    }

    for (unsigned int sub = 0; sub < l.number_of_segments(); ++sub) {
      for (int i = -1; i < D; ++i) {
        comp_lst[sub][i + 1].resize(l.number_of_runs(i, sub), -1);
        // total_number_of_runs+=l.number_of_runs(i,sub);
      }
    }

    bool is_first_run = true;
    int node_of_first_run = 0;
    int node_of_last_run = 0;

    // cycle through
    for (typename LStype::template const_iterator_neighbor_filtered<
             typename LStype::filter_all, 1>
             it(l);
         !it.is_finished(); it.next()) {

      int &tc = comp_lst[it.center().get_segment_num()][it.center().get_level()]
                        [it.center().run_type_position()];

      // -1 means it is not set yet
      if (tc == -1) {
        for (int k = 0; k < 2 * D; ++k) {
          const int &tn = comp_lst[it.neighbor(k).get_segment_num()]
                                  [it.neighbor(k).get_level()]
                                  [it.neighbor(k).run_type_position()];
          if (tn != -1) {
            // connected() { return it.sign() == it2.sign() }
            if (connected(it.center(), it.neighbor(k))) {
              tc = tn;
              break;
            }
          }
        }
      }

      // it is still not set, so add new vertex
      if (tc == -1) {
        tc = num_components;
        boost::add_vertex(Graph);
        ++num_components;
      }

      // check if edge can be set
      for (int k = 0; k < 2 * D; ++k) {
        int &tn = comp_lst[it.neighbor(k).get_segment_num()]
                          [it.neighbor(k).get_level()]
                          [it.neighbor(k).run_type_position()];
        if (connected(it.center(), it.neighbor(k))) {
          if (tn != -1) {
            if (tc != tn)
              boost::add_edge(tc, tn, Graph);
          } else {
            tn = tc;
          }
        }
      }

      // set special nodes
      if (is_first_run) {
        is_first_run = false;
        node_of_first_run = tc;
      }
      node_of_last_run = tc;
    }

    assert(boost::num_vertices(Graph) == num_components);
    std::vector<int> component(boost::num_vertices(Graph));

    unsigned int num_components_after =
        connected_components(Graph, &component[0]);

    // determine component number of source region

    int source_node = (is_open_boundary_negative) ? component[node_of_first_run]
                                                  : component[node_of_last_run];

    Connectivities.clear();
    for (typename LStype::template const_iterator_neighbor_filtered<
             typename LStype::filter_active, 1>
             it(l);
         !it.is_finished(); it.next()) {


      if (it.center().sign() == lvlset::POS_SIGN) {
        assert(it.center().get_segment_num() < l.number_of_segments());
        Connectivities.push_back(
            component[comp_lst[it.center().get_segment_num()][0]
                              [it.center().run_type_position()]] ==
            source_node); // TODO
      } else {
        int k;
        for (k = 0; k < 2 * D; ++k) {
          if (component[comp_lst[it.neighbor(k).get_segment_num()]
                                [it.neighbor(k).get_level()]
                                [it.neighbor(k).run_type_position()]] ==
              source_node)
            break;
        }
        Connectivities.push_back(k != 2 * D);
      }
    }

    return std::make_pair(num_components, num_components_after);
  }
};

#endif // LS_MARK_VOID_POINTS_HPP
