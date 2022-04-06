#include "basis.h"
vector<long> working_vertex;
vector<long> next_working_vertex;
vector<long> is_pending;
vector<long> vertex_to_removed;
void graph_reduction(int threshold) {
    for (auto v : remaining_vertex) {
        if (We[v] + vertex_neighbor_weight[v] <= threshold) {
            vertex_to_removed.push_back(v);
        }
    }
    while (!vertex_to_removed.empty()) {
        long i = *vertex_to_removed.rbegin();
        vertex_to_removed.pop_back();
        for (auto v:neighbor[i]) {
            vector<long>::size_type j = 0;
            for (; j<neighbor_len[v]; j++) {
                if (neighbor[v][j] == i) {
                    break;
                }
            }
            neighbor[v][j] = *neighbor[v].rbegin();
            neighbor[v].pop_back();
            vertex_neighbor_weight[v] -= We[i];
            neighbor_len[v]--;
			if (We[v] + vertex_neighbor_weight[v] + We[i] > threshold &&
					We[v] + vertex_neighbor_weight[v] <= threshold) {
				vertex_to_removed.push_back(v);
			}
        }
        //remove i
        neighbor[i].clear();
        remaining_vertex.remove(i);
    }
}

void graph_reduction_iterative(int threshold) {

    cout << "-------: graph reduction iterative" << endl;
	working_vertex.clear();  //clear
	next_working_vertex.clear();  //clear
	for (auto v : remaining_vertex) {
		working_vertex.push_back(v);
		is_pending[v] = true; //true if v in working_vertex or next_working_vertex
	}
    while (!working_vertex.empty()) {
		for (vector<long>::size_type i = 0; i < working_vertex.size(); ++i) {
			auto v = working_vertex[i];
			if (We[v] + vertex_neighbor_weight[v] <= threshold) {
				for (auto u : neighbor[v]) {
					vector<long>::size_type j = 0;
					for (; j < neighbor[u].size(); ++j) {
						if (neighbor[u][j] == v) {
							break;
						}
					}
					neighbor[u][j] = *neighbor[u].rbegin();
					neighbor[u].pop_back();
					vertex_neighbor_weight[u] -= We[v];
                    neighbor_len[u] --;
                    //Cy: vertex that need to be traversed again after traversing them. a vertex has been judged, some subsequent vertices adjacent to v are deleted, which in turn affects vertex v, so v needs to be added to next_working_vertex. is_pending[v] = false: denotes v has been judged.
                    //if u is not in working_vertex and next_working_vertex, add it to next_working_vertex
					if (!is_pending[u]) {
						next_working_vertex.push_back(u);
						is_pending[u] = true;
					}
				}
				neighbor[v].clear();
				remaining_vertex.remove(v);
			}
			is_pending[v] = false;
		}
		working_vertex.swap(next_working_vertex);
		next_working_vertex.clear();
	}
}
