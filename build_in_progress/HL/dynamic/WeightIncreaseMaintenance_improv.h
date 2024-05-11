#pragma once

#include <queue>
#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void SPREAD1(graph_v_of_v_idealID &instance_graph, vector<vector<two_hop_label_v1>> *L,
			 std::vector<affected_label> &al1, std::vector<pair_label> *al2, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{
	for (auto &al : al1)
	{
		int u = al.first;
		int v = al.second;
		double d_u = al.dis;
		std::queue<std::pair<int, double>> Q;
		Q.push({u, d_u});
		while (!Q.empty())
		{
			auto [x, d_x] = Q.front();
			Q.pop();
			auto temp = search_sorted_two_hop_label2((*L)[x], v);
			if(temp.second == -1) {
				insert_sorted_two_hop_label((*L)[x], v, MAX_VALUE);
				temp = search_sorted_two_hop_label2((*L)[x], v);
			}
			(*L)[x][temp.second].distance = MAX_VALUE;
			//(*L)[x][v].distance = MAX_VALUE;
			al2->push_back(pair_label(x, v));
			for (auto &x_n : instance_graph[x])
			{
				if (v < x_n.first)
				{
					double d_xn = d_x + x_n.second;
					if ((search_sorted_two_hop_label2((*L)[x_n.first], v).second != -1) && ((fabs(search_sorted_two_hop_label2((*L)[x_n.first], v).first - d_xn)) < 1e-9))
					{
						Q.push({x_n.first, d_xn});
					}
					else continue;
				}
			}
		}
	}
	/*TO DO 2*/
}

void SPREAD2(graph_v_of_v_idealID &instance_graph, vector<vector<two_hop_label_v1>> *L, PPR_type *PPR,
			 std::vector<pair_label> &al2, std::vector<affected_label> *al3, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{
	/*TO DO 3*/
	for (auto &it : al2)
	{ // it.first是x,it.second是y
		auto temp = PPR_retrieve(*PPR, it.first, it.second);
		
		//我需要遍历temp并上it.second，但不改变原先的PPR
		std::vector<int> temp2;
		for (auto t : temp)
		{
			temp2.push_back(t);
		}
		temp2.push_back(it.second);
		for (auto t : temp2)
		{
			if (t < it.first)
			{
				double dis = MAX_VALUE;
				for (auto &xn : instance_graph[it.first])
				{
					auto temp = search_sorted_two_hop_label2((*L)[xn.first], t);
					double new_dis;
					if(temp.second == -1) {
						continue;
					}
					else new_dis = temp.first + xn.second;
					dis = std::min(dis, new_dis);
				}
				if (dis < graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(*L, it.first, t))
				{
					al3->push_back(affected_label(it.first, t, dis));
				}
				else
				{
					auto [_, hub_v] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, it.first, t);
					PPR_insert(*PPR, it.first, hub_v, t);
					PPR_insert(*PPR, t, hub_v, it.first);
				}
			}
			if (it.first < t)
			{
				double dis = MAX_VALUE;
				for (auto &tn : instance_graph[t])
				{
					auto temp = search_sorted_two_hop_label2((*L)[tn.first], it.first);
					double new_dis;
					if(temp.second == -1) {
						continue;
					}
					else new_dis = temp.first + tn.second;
					dis = std::min(dis, new_dis);
				}
				if (graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(*L, t, it.first) > dis)
				{
					al3->push_back(affected_label(t, it.first, dis));
				}
				else
				{
					auto [_, hub_v] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, t, it.first);
					PPR_insert(*PPR, t, hub_v, it.first);
					PPR_insert(*PPR, it.first, hub_v, t);
				}
			}
		}
	}
}

void SPREAD3(graph_v_of_v_idealID &instance_graph, vector<vector<two_hop_label_v1>> *L, PPR_type *PPR, std::vector<affected_label> &al3,
			 ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{
	for (auto &al : al3)
	{
		int u = al.first;
		int v = al.second;
		double d_u = al.dis;
		double query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(*L, u, v);
		int common_hub = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, u, v).second;

		if (query_result <= d_u)
		{
			PPR_insert(*PPR, u, common_hub, v);
			PPR_insert(*PPR, v, common_hub, u);
			continue;
		}
		std::vector<double> Dis(instance_graph.size(), -1);
		Dis[u] = d_u;
		// Q为优先队列，存储(d_x, x)

		boost::heap::fibonacci_heap<std::pair<double, int>, boost::heap::compare<std::greater<std::pair<double, int>>>> fh;
		// std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> Q;
		std::unordered_map<int, boost::heap::fibonacci_heap<std::pair<double, int>, boost::heap::compare<std::greater<std::pair<double, int>>>>::handle_type> handles;
		// Q.push({d_u, u});
		auto handle = fh.push({d_u, u});
		handles[u] = handle;

		while (!fh.empty())
		{
			auto [d_x, x] = fh.top();
			fh.pop();

			//(*L)[x][v].distance = std::min(d_x, (*L)[x][v].distance);
			auto temp = search_sorted_two_hop_label2((*L)[x], v);
			if(temp.second == -1) {
				insert_sorted_two_hop_label((*L)[x], v, MAX_VALUE);
				temp = search_sorted_two_hop_label2((*L)[x], v);
			}
			(*L)[x][temp.second].distance = std::min(d_x, (*L)[x][temp.second].distance);

			for (auto &x_n : instance_graph[x])
			{
				if (v < x_n.first)
				{
					if (Dis[x_n.first] == -1)
					{
						Dis[x_n.first] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(*L, x_n.first, v);
					}
					double w_xn_x = x_n.second;
					if ( Dis[x_n.first] > d_x + w_xn_x)
					{
						if(w_xn_x != MAX_VALUE) {
							Dis[x_n.first] = d_x + w_xn_x;
						}
						else Dis[x_n.first] = MAX_VALUE;
						if (handles.find(x_n.first) != handles.end())
						{
							//{Dis[x_n.first], x_n.first}是一个斐波那契堆里的元素
							*(handles[x_n.first]) = std::make_pair(Dis[x_n.first], x_n.first);
							fh.update(handles[x_n.first]);
						}
						else
						{
							handles[x_n.first] = fh.push({Dis[x_n.first], x_n.first});
						}
					}
					else
					{
						common_hub = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, x_n.first, v).second;
						PPR_insert(*PPR, x_n.first, common_hub, v);
						PPR_insert(*PPR, v, common_hub, x_n.first);
					}
				}
			}
		}
	}

	
	/*TO DO 4*/
}

void WeightIncreaseMaintenance_improv(graph_v_of_v_idealID &instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1 &mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new,
									  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

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
