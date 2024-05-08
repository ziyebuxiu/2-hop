#pragma once

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void WeightDecreaseMaintenance_improv_step1(int v1, int v2, weightTYPE w_new, vector<vector<two_hop_label_v1>> *L, PPR_type *PPR, std::vector<affected_label> *CL,
											ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

	for (int sl = 0; sl < 2; sl++)
	{
		if (sl == 1)
		{
			swap(v1, v2);
		}
		for (auto it : (*L)[v1])
		{
			if (it.vertex <= v2)
			{
				results_dynamic.emplace_back(pool_dynamic.enqueue([it, v2, L, PPR, w_new, CL]
																  {

					auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, it.vertex, v2); // query_result is {distance, common hub}
					if (query_result.first > it.distance + w_new) {
						mtx_595_1.lock();
						CL->push_back(affected_label{ v2 , it.vertex, it.distance + w_new });
						mtx_595_1.unlock();
					}
					else {
						auto search_result = search_sorted_two_hop_label((*L)[v2], it.vertex);
						if (search_result > it.distance + w_new && search_result != MAX_VALUE) {
							mtx_595_1.lock();
							CL->push_back(affected_label{ v2, it.vertex, it.distance + w_new });
							mtx_595_1.unlock();
						}
						if (query_result.second != it.vertex) {
							mtx_5952[v2].lock();
							PPR_insert(*PPR, v2, query_result.second, it.vertex);
							mtx_5952[v2].unlock();
						}
						if (query_result.second != v2) {
							mtx_5952[it.vertex].lock();
							PPR_insert(*PPR, it.vertex, query_result.second, v2);
							mtx_5952[it.vertex].unlock();
						}
					}

					return 1; }));
			}
		}
	}

	for (auto &&result : results_dynamic)
	{
		result.get();
	}
	std::vector<std::future<int>>().swap(results_dynamic);
}

void DIFFUSE(graph_v_of_v_idealID &instance_graph, vector<vector<two_hop_label_v1>> *L, PPR_type *PPR, std::vector<affected_label> &CL,
			 ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

	/*TO DO 1*/
	// for (auto &cl : CL)
	// {
	// 	int u = cl.first, v = cl.second;
	// 	weightTYPE du = cl.dis;

	// 	auto task = [u, v, du, &instance_graph, L, PPR]()
	// 	{
	// 		std::vector<int> Dis(instance_graph.size(), -1);
	// 		std::queue<std::pair<int, weightTYPE>> Q;
	// 		Dis[u] = du;
	// 		Q.push({u, du});
	// 		int common_hub = -1;
	// 		while (!Q.empty())
	// 		{
	// 			auto [x, dx] = Q.front();
	// 			Q.pop();

	// 			// Updating the two-hop label for node x
	// 			if ((*L)[x][v].distance > dx)
	// 			{
	// 				(*L)[x][v].distance = dx;
	// 			}

	// 			// Exploring neighbors of x
	// 			for (auto &xn : instance_graph[x])
	// 			{
	// 				if (xn.first > v)
	// 				{
	// 					int neighbor = xn.first;
	// 					weightTYPE weight_x_xn = xn.second;
	// 					weightTYPE new_distance = dx + weight_x_xn;

	// 					if (Dis[neighbor] == -1)
	// 					{
	// 						Dis[neighbor] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(*L, neighbor, v);
	// 					}
	// 					if(Dis[neighbor] > new_distance){
	// 						Dis[neighbor] = new_distance;
	// 						Q.push({neighbor, new_distance});
	// 					}
	// 					// Update PPR if applicable
	// 					else if (v != neighbor && std::min((*L)[neighbor][v].distance, new_distance) > dx + weight_x_xn)
	// 					{
	// 						mtx_5952[neighbor].lock();
	// 						PPR_insert(*PPR, neighbor, v, new_distance);
	// 						mtx_5952[neighbor].unlock();
	// 					}
	// 					// (*PPR)[neighbor][]
	// 				}
	// 				else
	// 					continue;
	// 			}
	// 		}
	// 		return 1;
	// 	};

	// 	results_dynamic.emplace_back(pool_dynamic.enqueue(task));
	// }

	// // Wait for all tasks to complete
	// for (auto &&result : results_dynamic)
	// {
	// 	result.get();
	// }
	// results_dynamic.clear(); // Clear the futures to be ready for new tasks if needed
}

void WeightDecreaseMaintenance_improv(graph_v_of_v_idealID &instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1 &mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new,
									  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

	std::vector<affected_label> CL;
	WeightDecreaseMaintenance_improv_step1(v1, v2, w_new, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic);

	DIFFUSE(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic);
}
