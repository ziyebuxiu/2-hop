#pragma once

#include <queue>
#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void SPREAD1(graph_v_of_v_idealID &instance_graph, vector<vector<two_hop_label_v1>> *L,
			 std::vector<affected_label> &al1, std::vector<pair_label> *al2, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

	// /*TO DO 2*/
	for (const auto &item : al1)
	{
		int u = item.first, v = item.second;
		weightTYPE du = item.dis;
		results_dynamic.emplace_back(pool_dynamic.enqueue([u, v, du, L, al2, &instance_graph]()
														  {
            std::queue<std::pair<int, weightTYPE>> Queue;
            Queue.push({u, du});
            while (!Queue.empty()) {
                auto [x, dx] = Queue.front();
                Queue.pop();
                (*L)[x][v].distance = std::numeric_limits<weightTYPE>::infinity();  // Invalidate this label
                al2->push_back(pair_label(x, v));  // Prepare for further spreading
                for (auto &xn : instance_graph[x]) {
					if(v > xn.first){
						if (((*L)[xn.first][v].distance == instance_graph[x][xn.first].second + dx)) {
                        	Queue.push({xn.first, dx + instance_graph[x][xn.first].second});
                    	}
					}

                }
            }
            return 1; }));
	}
}

void SPREAD2(graph_v_of_v_idealID &instance_graph, vector<vector<two_hop_label_v1>> *L, PPR_type *PPR,
			 std::vector<pair_label> &al2, std::vector<affected_label> *al3, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

	// /*TO DO 3*/
	for (const auto &pair : al2)
	{
		auto task = [&instance_graph, L, PPR, al3, pair]()
		{
			int x = pair.first, y = pair.second;
			for (auto &target : PPR[x][y])
			{
				// target.first和second分别代表什么?
				weightTYPE min_distance = std::numeric_limits<weightTYPE>::infinity();
				// int common_hub = -1;
				if (target.first > x)
				{
					// Calculate the potential new shortest path through neighbors
					for (auto &neighbor : instance_graph[x])
					{
						auto [dnew, index] = search_sorted_two_hop_label2((*L)[neighbor.first], target.first);
						if(index == -1) continue;
						weightTYPE new_distance = dnew + instance_graph[x][neighbor.first].second;
						if (new_distance < min_distance)
						{
							min_distance = new_distance;
							// common_hub = neighbor.first;
						}
					}
					auto [dis, hc] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, x, target.first);
					if (min_distance < dis)
					{
						al3->push_back(affected_label(x, target.first, min_distance));
					}
					else
					{
						(*PPR)[x][hc].second.push_back(target.first);
						(*PPR)[target.first][hc].second.push_back(x);
					}
				}
				if (x > target.first)
				{
					// Calculate the potential new shortest path through neighbors
					for (auto &neighbor : target.second)
					{
						auto [dnew, index] = search_sorted_two_hop_label2((*L)[neighbor], x);
						if(index == -1) continue;
						weightTYPE new_distance = dnew + instance_graph[target.first][neighbor].second;
						// weightTYPE new_distance = (*L)[neighbor][x].distance + instance_graph[target.first][neighbor].second;
						if (new_distance < min_distance)
						{
							min_distance = new_distance;
							// common_hub = neighbor;
						}
					}
					auto [dis, hc] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, x, target.first);
					if (min_distance < dis)
					{
						al3->push_back(affected_label(target.first, x, min_distance));
					}
					else
					{
						(*PPR)[target.first][hc].second.push_back(x);
						(*PPR)[x][hc].second.push_back(target.first);
					}
				}
			}
			return 1;
		};
		results_dynamic.emplace_back(pool_dynamic.enqueue(task));
	}
}

void SPREAD3(graph_v_of_v_idealID &instance_graph, vector<vector<two_hop_label_v1>> *L, PPR_type *PPR, std::vector<affected_label> &al3,
			 ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

	/*TO DO 4*/
	// for (auto &affected : al3) {
	//     int u = affected.first, v = affected.second;
	//     weightTYPE du = affected.dis;
	//     results_dynamic.emplace_back(pool_dynamic.enqueue([u, v, du, L, PPR, &instance_graph]() {
	//         std::queue<std::pair<int, weightTYPE>> Q;
	//         Q.push({u, du});
	//         std::vector<weightTYPE> Dis(instance_graph.size(), std::numeric_limits<weightTYPE>::max());
	//         Dis[u] = du;

	//         while (!Q.empty()) {
	//             auto [x, dx] = Q.front();
	//             Q.pop();
	//             (*L)[x][v].distance = std::min((*L)[x][v].distance, dx);  // Update the label
	//             for (auto &neighbor : instance_graph[x]) {
	//                 if (neighbor.first > v && Dis[neighbor.first] > dx + neighbor.second) {
	//                     Dis[neighbor.first] = dx + neighbor.second;
	//                     Q.push({neighbor.first, Dis[neighbor.first]});
	//                     // Update PPR if the condition matches
	//                     PPR->at(neighbor.first)[v].second.push_back(x);
	//                     PPR->at(x)[v].second.push_back(neighbor.first);
	//                 }
	//             }
	//         }
	//         return 1;
	//     }));
	// }
}

void WeightIncreaseMaintenance_improv(graph_v_of_v_idealID &instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1 &mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new,
									  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{
	// 这个函数用于处理边权增加的情况，即w_old < w_new，实现的功能是更新两个节点之间的所有两跳标签
	std::vector<affected_label> al1, al3;
	std::vector<pair_label> al2;

	/*it's slow to paralize the following part*/
	for (auto it : mm.L[v1])
	{
		if (it.vertex <= v2 && abs(search_sorted_two_hop_label(mm.L[v2], it.vertex) - it.distance - w_old) < 1e-5)
		{
			al1.push_back(affected_label(v2, it.vertex, it.distance + w_old));
		}
	}
	for (auto it : mm.L[v2])
	{
		if (it.vertex <= v1 && abs(search_sorted_two_hop_label(mm.L[v1], it.vertex) - it.distance - w_old) < 1e-5)
		{
			al1.push_back(affected_label(v1, it.vertex, it.distance + w_old));
		}
	}

	// cout << "al1.size() " << al1.size() << endl;

	SPREAD1(instance_graph, &mm.L, al1, &al2, pool_dynamic, results_dynamic);
	SPREAD2(instance_graph, &mm.L, &mm.PPR, al2, &al3, pool_dynamic, results_dynamic);
	SPREAD3(instance_graph, &mm.L, &mm.PPR, al3, pool_dynamic, results_dynamic);

	// for (auto it : al2) {
	//	cout << "al2 " << it.first << " " << it.second << endl;
	// }
	// for (auto it : al3) {
	//	cout << "al3 " << it.first << " " << it.second << " " << it.dis << endl;
	// }
}
