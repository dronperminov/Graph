#include <iostream>
#include <iomanip>

#include "Graph.hpp"

using namespace std;

int main() {
	Graph graph(7, false);

	graph.AddEdge(0, 3, 1);
	graph.AddEdge(0, 4, 2);
	graph.AddEdge(3, 5, 3);
	graph.AddEdge(3, 6, 4);
	graph.AddEdge(4, 6, 5);
	graph.AddEdge(5, 6, 6);
	graph.AddEdge(1, 4, 7);
	graph.AddEdge(1, 2, 8);
	graph.AddEdge(1, 5, 9);

	cout << "Matrix of adjacency: " << endl;
	graph.PrintAdjacencyMatrix();

	cout << endl << "Matrix of incidence: " << endl;
	graph.PrintIncidenceMatrix();

	cout << endl << "List of edges: " << endl;
	graph.PrintEdgesList();

	cout << endl << "List of connectivities: " << endl;
	graph.PrintConnectivityList();

	cout << endl << "DFS from 0: ";
	graph.DFS(0);

	cout << endl << "DFS from 2: ";
	graph.DFS(2);
	cout << endl;

	cout << endl << "BFS from 0: ";
	graph.BFS(0);

	cout << endl << "BFS from 2: ";
	graph.BFS(2);
	cout << endl;

	if (graph.IsConnected()) {
		cout << endl << "Graph is connected" << endl;
	}
	else {
		cout << endl << "Graph is not connected" << endl;
	}

	if (graph.IsVertexReachable(0, 2)) {
		cout << endl << "Vertex 0 is reachable from vertex 2" << endl;
	}
	else {
		cout << endl << "Vertex 0 is not reachable from vertex 2" << endl;	
	}

	/************************* DIJKSTRA *************************/
	cout << endl << "Try to find way from 1 to 3:" << endl;
	std::vector<int> way = graph.DijkstraAlgorithm(1, 3);

	if (way.size() == 0) {
		cout << "No way from 1 to 3" << endl;
	}
	else {
		int distance = 0;
		// выводим кратчайщий путь
		cout << "Way:";

		for (int i = way.size() - 1; i >= 0; i--) {
			if (i > 0)
				distance += graph.GetWeight(way[i], way[i - 1]);

			cout << " " << way[i];
		}

		cout << endl;
		cout << "Distance: " << distance << endl;
	}

	/************************* MST *************************/
	cout << endl << "MST: ";
	vector<int> tree = graph.MST();

	for (size_t i = 0; i < tree.size(); i += 2)
		cout << tree[i] << "-" << tree[i + 1] << " ";
	cout << endl;

	/************************* Floid-Yorshall *************************/
	vector<vector<int>> W = graph.FloidYorshallAlgorithm();

	cout << endl << "Matrix of Floid-Yorshall algorithm:" << endl;
	for (size_t i = 0; i < W.size(); i++) {
		for (size_t j = 0; j < W[i].size(); j++) {
			if (W[i][j] == INF) {
				cout << "  -  ";
			}
			else {
				cout << setw(3) << W[i][j] << " ";
			}
		}

		cout << endl;
	}

	/************************* Vertices degrees *************************/
	cout << endl << "Degrees of vertices:" << endl;
	for (int i = 0; i < 7; i++) {
		cout << "deg of v" << i << ": " << graph.GetVertexDegree(i) << endl;
	}

	/************************* Reachability matrix *************************/
	vector<vector<bool>> matrix = graph.GetReachabilityMatrix();
	cout << endl << "Reachability matrix: " << endl;

	for (size_t i = 0; i < matrix.size(); i++) {
		for (size_t j = 0; j < matrix[i].size(); j++)
			cout << setw(3) << matrix[i][j] << " ";

		cout << endl;
	}

}