#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <queue>
#include <vector>

const int INF = 10000000;

// класс графа, внутреннее представление - матрица смежности
class Graph {
	// стрктура ребра
	struct Edge {
		int weight;
		int v1;
		int v2;
	};

	int vertices; // число вершин
	bool isDirected; // ориентированный ли граф
	int **matrix; // матрица смежности

	void PrintDFS(int u, std::vector<bool> &visited) const; // рекурсивная функция вывода графа при обходе в глубину
	int VisitedVerticesCount(int u, std::vector<bool> &visited) const; // количество посещённых вершин из вершины u
	void DFS(int u, std::vector<bool> &visited) const; // рекурсивная функция обхода графа в глубину 

public:
	Graph(int vertices, bool isDirected = true); // создание графа из числа вершин

	void AddEdge(int v1, int v2, int weight = 1); // добавление ребра v1-v2
	bool HaveEdge(int v1, int v2) const; // проверка наличия ребра v1-v2
	int GetWeight(int v1, int v2) const; // получение веса ребра v1-v2

	void PrintAdjacencyMatrix() const; // вывод матрицы смежности
	void PrintIncidenceMatrix() const; // вывод матрицы индидентности
	void PrintEdgesList() const; // вывод списка рёбер
	void PrintConnectivityList() const; // вывод списка связности

	void DFS(int u) const; // обход в глубину из вершины u
	void BFS(int u) const; // обход в ширину из вершины u

	int GetVertexDegree(int v) const; // получение степени вершины

	bool IsConnected() const; // является ли граф связным
	bool IsVertexReachable(int v, int u) const; // проверка достижимости вершины v из вершины u

	std::vector<int> DijkstraAlgorithm(int u, int v) const; // поиск кратчайшего пути, алгоритм Дейкстры
	std::vector<std::vector<int>> FloidYorshallAlgorithm() const; // алгоритм Флойда-Уоршала
	std::vector<std::vector<bool>> GetReachabilityMatrix() const; // получение матрицы достижимости
	std::vector<int> PrimMST() const; // построение минимального остовного дерева, алгоритм Прима
	std::vector<int> KruskalMST() const; // построение минимального остовного дерева, алгоритм Крускала

	~Graph(); // деструктор
};

Graph::Graph(int vertices, bool isDirected) {
	this->vertices = vertices;
	this->isDirected = isDirected;

	matrix = new int*[vertices]; // выделяем память под матрцу

	for (int i = 0; i < vertices; i++) {
		matrix[i] = new int[vertices];

		// обнуляем матрицу (изначально нет рёбер)
		for (int j = 0; j < vertices; j++)
			matrix[i][j] = 0;
	}
}

// рекурсивная функция вывода графа при обходе в глубину
void Graph::PrintDFS(int u, std::vector<bool> &visited) const {
	visited[u] = true; // помечаем вершину как посещённую

	std::cout << " " << u; // выводим посещаемую вершину

	// проходимся по смежным с u вершинам
	for (int v = 0; v < vertices; v++) {
		if (matrix[u][v] == 0)
			continue;

		// если вершина не была посещена
		if (!visited[v]) {
			PrintDFS(v, visited); // посещаем её
		}
	}
}

// обход графа в глубину
void Graph::DFS(int u, std::vector<bool> &visited) const {
	visited[u] = true; // помечаем вершину как посещённую

	// проходимся по смежным с u вершинам
	for (int v = 0; v < vertices; v++) {
		if (matrix[u][v] == 0)
			continue;

		// если вершина не была посещена
		if (!visited[v]) {
			DFS(v, visited); // посещаем её
		}
	}
}

// количество посещённых вершин из вершины u
int Graph::VisitedVerticesCount(int u, std::vector<bool> &visited) const {
	int visitedVertices = 1;
	
	visited[u] = true; // помечаем вершину как посещённую

	// проходимся по смежным с u вершинам
	for (int v = 0; v < vertices; v++) {
		if (matrix[v][u] == 0)
			continue;

		// если вершина не была посещена
		if (!visited[v])
			visitedVertices += VisitedVerticesCount(v, visited); // посещаем её
	}

	return visitedVertices; // возвращаем количество посещённых вершин
}

// добавление ребра
void Graph::AddEdge(int v1, int v2, int weight) {
	// проверяем на корректность значения вершин для ребра
	if (v1 < 0 || v2 < 0 || v1 >= vertices || v2 >= vertices)
		throw std::string("Graph::AddEdge - incorrect vertex(es)");

	matrix[v1][v2] = weight; // записываем ребро в матрицу

	if (!isDirected)
		matrix[v2][v1] = weight; // добавляем ребро v2-v1, если граф неориентированный
}

// проверка наличия ребра v1-v2
bool Graph::HaveEdge(int v1, int v2) const {
	// проверяем на корректность значения вершин для ребра
	if (v1 < 0 || v2 < 0 || v1 >= vertices || v2 >= vertices)
		throw std::string("Graph::HaveEdge - incorrect vertex(es)");

	return matrix[v1][v2] != 0; // ребро есть, если в матрице не ноль
}

// получение веса ребра v1-v2
int Graph::GetWeight(int v1, int v2) const {
	// проверяем на корректность значения вершин для ребра
	if (v1 < 0 || v2 < 0 || v1 >= vertices || v2 >= vertices)
		throw std::string("Graph::HaveEdge - incorrect vertex(es)");

	return matrix[v1][v2];
}

// вывод матрицы смежности
void Graph::PrintAdjacencyMatrix() const {
	for (int i = 0; i < vertices; i++) {
		for (int j = 0; j < vertices; j++)
			std::cout << std::setw(3) << matrix[i][j] << " ";
		
		std::cout << std::endl;
	}
}

// вывод матрицы индидентности
void Graph::PrintIncidenceMatrix() const {
	int edges = 0;

	// считаем число рёбер
	for (int i = 0; i < vertices; i++)
		for (int j = isDirected ? 0 : i; j < vertices; j++)
			if (matrix[i][j] != 0)
				edges++;

	// выделяем память под матрицу инцидентности
	int **M = new int*[vertices];

	for (int i = 0; i < vertices; i++) {
		M[i] = new int[edges];

		for (int j = 0; j < edges; j++)
			M[i][j] = 0;
	}

	int index = 0;

	// заполняем матрицу инцидентности
	for (int i = 0; i < vertices; i++) {
		for (int j = isDirected ? 0 : i; j < vertices; j++) {
			if (matrix[i][j] != 0) {
				M[i][index] = 1;
				M[j][index] = 1;
				index++;
			}
		}
	}

	// выводим матрицу инцидентности
	for (int i = 0; i < vertices; i++) {
		for (int j = 0; j < edges; j++)
			std::cout << " " << M[i][j];
		
		std::cout << std::endl;
		delete[] M[i];
	}

	delete[] M;
}

// вывод списка рёбер
void Graph::PrintEdgesList() const {
	for (int i = 0; i < vertices; i++)
		for (int j = isDirected ? 0 : i; j < vertices; j++)
			if (matrix[i][j] != 0)
				std::cout << i << " " << j << std::endl;
}

// вывод списка связности
void Graph::PrintConnectivityList() const {
	for (int i = 0; i < vertices; i++) {
		std::cout << i << ":";

		for (int j = 0; j < vertices; j++)
			if (matrix[i][j] != 0)
				std::cout << " " << j;

		std::cout << std::endl;
	}
}

void Graph::DFS(int u) const {
	std::vector<bool> visited(vertices, false); // вектор флагов посещённости

	PrintDFS(u, visited); // выполняем рекурсивную функцию обхода в глубину
}

void Graph::BFS(int u) const {
	std::queue<int> q;
	std::vector<bool> visited(vertices, false); // массив флагов посещения вершин

	q.push(u);
	visited[u] = true;

	std::cout << " " << u;

	while (!q.empty()) {
		int v = q.front();
		q.pop();

		for (int i = 0; i < vertices; i++) {
			if (matrix[v][i] != 0 && !visited[i]) {
				visited[i] = true;
				q.push(i);

				std::cout << " " << i;
			}
		}
	}
}

// получение степени вершины
int Graph::GetVertexDegree(int v) const {
	if (v < 0 || v >= vertices)
		throw std::string("Graph::GetVertexDegree() - incorrect vertex");

	int degree = 0;

	for (int i = 0; i < vertices; i++)
		if (matrix[v][i] != 0)
			degree++;

	return degree;
}

// является ли граф связным
bool Graph::IsConnected() const {
	std::vector<bool> visited(vertices, false); // вектор флагов посещённости

	return VisitedVerticesCount(0, visited) == vertices; // граф связный, если все вершин могут быть посещены
}

// проверка достижимости вершины v из вершины u
bool Graph::IsVertexReachable(int v, int u) const {
	std::vector<bool> visited(vertices, false); // вектор флагов посещённости

	DFS(u, visited); // выполняем обход в глубину из вершины u

	return visited[v];
}

// поиск кратчайшего пути, алгоритм Дейкстры
std::vector<int> Graph::DijkstraAlgorithm(int s, int p) const {
	std::vector<bool> wasFound(vertices, false); // вектор флагов, была ли посещена конкретная точка
	std::vector<int> distance(vertices, INF); // вектор расстояний от точки s до любой другой

	distance[s] = 0; // расстояние от s до s равно нулю
	wasFound[s] = true; // стартовая вершина была посещена

	int last = s; // последняя посещённая вершина - s

	for (int i = 0; i < vertices; i++) {
		for (int j = 0; j < vertices; j++) {
			if (matrix[last][j] != 0 && !wasFound[j]) {
				if (distance[j] > distance[last] + matrix[last][j])
					distance[j] = distance[last] + matrix[last][j];
			}
		}

		int minDist = INF; // считаем, что минимальное расстояние - бесконечность

		for (int j = 0; j < vertices; j++) {
			if (!wasFound[j] && minDist > distance[j]) {
				minDist = distance[j];
				last = j;
			}
		}

		wasFound[last] = true; // отмечаем вершину как посещённую
	}

	std::vector<int> way; // вектор пути

	// если полученноое расстояние равно бесконечности, значит пути нет
	if (distance[p] == INF)
		return way;

	// восстанавливаем путь
	int end = p; // конечная вершина

	way.push_back(end);
	int weight = distance[end];

	// пока не придём в начало
	while (end != s) {
		for (int i = 0; i < vertices; i++) { // просматриваем все вершины
			if (matrix[i][end] != 0) {  // если связь есть
				int tmp = weight - matrix[i][end]; // определяем вес пути из предыдущей вершины

				if (tmp == distance[i]) { // если вес совпал с рассчитанным, значит из этой вершины и был переход
					weight = tmp; // сохраняем новый вес
					end = i;       // сохраняем предыдущую вершину
					way.push_back(i); // и записываем ее в массив
				}
			}
		}
	}

	return way;
}

// алгоритм Флойда-Уоршала
std::vector<std::vector<int>> Graph::FloidYorshallAlgorithm() const {
	std::vector<std::vector<int>> W(vertices, std::vector<int>(vertices));

	for (int i = 0; i < vertices; i++) {
		for (int j = 0; j < vertices; j++) {
			if (matrix[i][j] == 0 && i != j) { // если связи нет
				W[i][j] = INF; // то расстояние бесконечно
			}
			else {
				W[i][j] = matrix[i][j]; // иначе копируем значение расстояния
			}
		}
	}

	// выполняем алгоритм Флойда-Уоршелла
	for (int k = 0; k < vertices; k++) {
		for (int i = 0; i < vertices; i++) {
			for (int j = 0; j < vertices; j++) {
				W[i][j] = W[i][j] < W[i][k] + W[k][j] ? W[i][j] : W[i][k] + W[k][j];
			}
		}
	}

	return W; // возвращаем матрицу
}

// получение матрицы достижимости
std::vector<std::vector<bool>> Graph::GetReachabilityMatrix() const {
	std::vector<std::vector<bool>> M(vertices, std::vector<bool>(vertices));

	for (int i = 0; i < vertices; i++)
		for (int j = 0; j < vertices; j++)
			M[i][j] = matrix[i][j] != 0;

	for (int k = 0; k < vertices; k++)
		for (int i = 0; i < vertices; i++)
			for (int j = 0; j < vertices; j++)
				M[i][j] = M[i][j] || (M[i][k] && M[k][j]);

	return M;
}

// построение минимального остовного дерева, алгоритм Прима
std::vector<int> Graph::PrimMST() const {
	if (isDirected)
		throw std::string("Graph::PrimMST() - graph is directed");

	std::vector<bool> used(vertices, false); // массив флагов посещения вершин
	std::vector<int> min_edge(vertices, INF); // массив минимальных длин рёбер
	std::vector<int> end_edge(vertices, -1); // массив концов вершин в остовном дереве
	std::vector<int> tree;

	min_edge[0] = 0; // начинаем с нулевой вершины

	for (int i = 0; i < vertices; i++) {
		int v = -1;

		// ищем вершину с минимальным весом, которая ещё не была посещена
		for (int j = 0; j < vertices; j++)
			if (!used[j] && (v == -1 || min_edge[j] < min_edge[v]))
				v = j;
		
		used[v] = true; // помечаем её как посещённую

		if (end_edge[v] != -1) {
			tree.push_back(v);
			tree.push_back(end_edge[v]);
		}

		for (int to = 0; to < vertices; to++) {
			if (matrix[v][to] != 0 && matrix[v][to] < min_edge[to]) {
				min_edge[to] = matrix[v][to];
				end_edge[to] = v;
			}
		}
	}

	return tree; // возвращаем дерево
}

// построение минимального остовного дерева, алгоритм Крускала
std::vector<int> Graph::KruskalMST() const {
	if (isDirected)
		throw std::string("Graph::KruskalMST() - graph is directed");

	std::vector<Edge> edges; // вектор рёбер

	// формируем вектор рёбер
	for (int i = 0; i < vertices; i++)
		for (int j = i; j < vertices; j++)
			if (matrix[i][j] != 0)
				edges.push_back({ matrix[i][j], i, j });

	// сортируем рёбра по возрастанию весов
	for (int k = edges.size() / 2; k > 0; k /= 2) {
		for (size_t i = k; i < edges.size(); i++) {
			int j = i;
			Edge tmp = edges[i];

			while (j >= k && tmp.weight < edges[j - k].weight) {
				edges[j] = edges[j - k];
				j -= k;
			}

			edges[j] = tmp;
		}
	}

	std::vector<int> treeId(vertices);

	for (int i = 0; i < vertices; i++)
		treeId[i] = i;

	std::vector<int> tree; // вектор для минимального остовного дерева

	for (size_t i = 0; i < edges.size(); i++) {
		int v1 = edges[i].v1;
		int v2 = edges[i].v2;

		if (treeId[v1] != treeId[v2]) {
			tree.push_back(v1);
			tree.push_back(v2);

			int oldId = treeId[v2];
			int newId = treeId[v1];

			for (int j = 0; j < vertices; j++)
				if (treeId[j] == oldId)
					treeId[j] = newId;
		}
	}

	return tree; // возвращаем дерево
}

// деструктор
Graph::~Graph() {
	// освобождаем память из-под матрицы
	for (int i = 0; i < vertices; i++)
		delete[] matrix[i];

	delete[] matrix;
}