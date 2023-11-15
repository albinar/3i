#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <cfloat>

#include <Rinternals.h>
#include <R.h>

#include "graehl/graph.h"
#include "numberGenerator.h"
#include "exceptions.h"
#include "dataRead.h"
#include "mystring.h"

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

using namespace std;
using namespace tfl;

int nonZeroCount(string input)
{
	int count = 0;
	for (int i = 0; input.c_str()[i]; i++)
		if (input.c_str()[i] != '0')
			count++;
	return count;
}

// adds to graph like dfs.
void addToGraph(Graph &g, map<string, int> &graphIdMap, map<string, double> &pvals, NumberGenerator &n, int &id, vector<bool> &visited)
{
//	cerr << "called." << endl;
	if (graphIdMap.find(n.text()) == graphIdMap.end())
	{
		graphIdMap[n.text()] = id;
		//Rprintf("befor loop, num: %s\tid: %d\n", n.text(), id);
		visited.push_back(false);
		id++;
	}

	for (int i = 0; i < n.nonZeroCount(); i++)
	{
		GraphArc a;
		NumberGenerator neighbor = n.neighbor(i);

		if (graphIdMap.find(neighbor.text()) == graphIdMap.end())
		{
			graphIdMap[neighbor.text()] = id;
			//Rprintf("loop, num: %s\tid: %d\n", neighbor.text(), id);
			visited.push_back(false);
			id++;
		}

		a.dest = graphIdMap[n.text()];
		a.source = graphIdMap[neighbor.text()];
		a.weight = -pvals[n.text()];
		a.data = NULL;
		g.states[a.source].arcs.push(a);
//		printf("(%s %s %lg)\n", neighbor.text(), n.text(), a.weight);
//		cerr << "source: " << n.text() << " dest: " << neighbor.text();
//		cerr << " one item added: " << a << endl;
//		cerr << g;

		if (!visited[graphIdMap[neighbor.text()]])
		{
//			fprintf(stderr, "----calling in this state: %s to %s\n", n.text(), neighbor.text());
			addToGraph(g, graphIdMap, pvals, neighbor, id, visited);
		}
	}
	visited[graphIdMap[n.text()]] = true;
}

// two first inputs are to be filled and should be empty as input
void createGraph(Graph &g, map<string, int> &graphIdMap, map<string, double> &pvals, NumberGenerator dest)
{
	g.nStates = pow(2.0, dest.nonZeroCount());
    g.states = new GraphState[g.nStates];
	
	graphIdMap.clear();
	int id = 0;

	vector<bool> visited;

	addToGraph(g, graphIdMap, pvals, dest, id, visited);
}

string getKey(map<string, int> m, int v)
{
	for (map<string, int>::iterator it = m.begin(); it != m.end(); it++)
		if (it->second == v)
			return it->first;
	return "NULL";
}

void addEdge(map<string, map<string, double> > &edge_list, string source, string dest, double value)
{
	if (edge_list.find(source) == edge_list.end())
		edge_list[source] = map<string, double>();
	if (edge_list[source].find(dest) == edge_list[source].end())
		edge_list[source][dest] = value;
}

void addNodeSize(map<string, double> &node_size, string node, double size)
{
	if (node_size.find(node) == node_size.end())
		node_size[node] = size;
	else
		node_size[node] += size;
}

string getHex(unsigned char value)
{
	char c1 = value / 16;
	char c2 = value % 16;
	if (c1 < 10)
		c1 = '0' + c1;
	else
		c1 = 'A' + c1 - 10;

	if (c2 < 10)
		c2 = '0' + c2;
	else
		c2 = 'A' + c2 - 10;
	string res;
	res.append(1, c1);
	res.append(1, c2);
	return res;
}

string diffPhenotype(const string &s1, const string &s2)
{
	string res = "";
	for (unsigned int i = 0; i < s1.size(); i++)
		if (s1[i] != s2[i])
		{
			if (s1[i] != '0')
				res += s1[i];
			else
				res += s2[i];
		}
		else
			res += '0';
	return res;
}

int get_trim_level(map<string, double> &pvals, map<string, int> &graphIdMap, List<GraphArc *> &p, int do_dynamic_trim, int dynamic_trim_tolerance, int do_static_trim, int static_trim_level)
{
  //Rprintf("__\n");
	int level=1;
	int increase_start_level = 0;
	int increasing = 0;
	double increase_start_score = 0;
	GraphArc *a;
	for ( ListIter<GraphArc *> arcIter(p) ; arcIter ; ++arcIter )
	{
		a = arcIter.data();
		if (level >= static_trim_level && do_static_trim)
		{
			if (increasing && do_dynamic_trim)
				return increase_start_level;
			else
				return level;
		}

		if (pvals[getKey(graphIdMap, a->dest)] < pvals[getKey(graphIdMap, a->source)] && !increasing)
		{
		  /*Rprintf("dest %lg\tsource %lg\tincreasing %d\tlevel %d\t%s\n", 
					pvals[getKey(graphIdMap, a->dest)],
					pvals[getKey(graphIdMap, a->source)],
					increasing,
					level,
					getKey(graphIdMap, a->source).c_str());*/
			increase_start_score = pvals[getKey(graphIdMap, a->source)];
			increase_start_level = level;
			increasing = 1;
		}

		if(increasing && pvals[getKey(graphIdMap, a->source)] > increase_start_score)
		{
			increasing = 0;
		}

		if (increasing && level - increase_start_level >= dynamic_trim_tolerance && do_dynamic_trim)
			return increase_start_level;

		level++;
	}
	return level;
}

extern "C"
{
	void c_analyze(char **s_signs, 
			int *signs_row_num, 
			double *s_pvals, 
			double *s_max_score,
			double *s_min_score,
			char **s_phenotype_set, 
			int *s_best_count, 
			int *trim_paths,
			int *trim_tolerance,
			int *trim_level,
			//outputs
			char **best_paths, 
			double *node_pvalue,
			char **gnodes, 
			char **gedges,
			int *node_count,
			int *edge_count)
	{
		map<string, double> pvals;
		map<string, string> signs;

		double max_score = *s_max_score;
		double min_score = *s_min_score;
		int data_count = *signs_row_num;
		int best_count = *s_best_count;
		int static_trim_level = *trim_level;
		int dynamic_trim_tolerance = *trim_tolerance;
		int do_static_trim = 0;
		if (static_trim_level > 0)
			do_static_trim = 1;
		int do_dynamic_trim = *trim_paths;

		//printf("****************\ndo trim: %d\n*****************\n", do_trim);
		string phenotype_set = *s_phenotype_set;

		for (int i = 0; i < data_count; i++)
		{
			//string sign_name = s_signs[i];
			double pval = s_pvals[i];
			string key = s_signs[i];
			//for (int j = 0; j < protein_count; j++)
			//{
			//	char tmp = '0' + s_values[i * protein_count + j];
			//	key += tmp;
			//}
			//printf("%s\t%s\t%lg\t\t", key.c_str(), sign_name.c_str(), pval);
			//inversZeroOne(key);
			pvals[key] = pval;
			//signs[key] = sign_name;
			//printf("%s\t%s\t%lg\n", key.c_str(), sign_name.c_str(), pval);
		}
		
		char *back_trace_start = *s_phenotype_set;
		NumberGenerator dest(strlen(back_trace_start), 3, back_trace_start);
		string start;
		start.insert(0, strlen(back_trace_start), '0');
		NumberGenerator source(strlen(back_trace_start), 3, start.c_str());
		Graph graph;
		map<string, int> graphIdMap;
		createGraph(graph, graphIdMap, pvals, dest);
		//Rprintf("%s", graph2str(graph).c_str());
		map<string, double> node_size;
		map<string, map<string, double> > edge_list;
		node_size.clear();
		edge_list.clear();
		//	cout << source.text() << " " << graphIdMap[source.text()] << "\n" << dest.text() << " " << graphIdMap[dest.text()] << "\n";
		//	return;

		int node_pos = 0;
		double min_edge = DBL_MAX;
		double max_edge = DBL_MIN;
		//printf("source:%s\ndest:%s\n", source.text(), dest.text());
		List<List<GraphArc *> > *paths = bestPaths(graph, graphIdMap[source.text()], graphIdMap[dest.text()], best_count);
		for ( ListIter<List<GraphArc *> > pathIter(*paths) ; pathIter ; ++pathIter ) 
		{
			double pathWeight = 0;
			GraphArc *a;
			int level = 1;
			int path_trim_level = get_trim_level(pvals, graphIdMap, pathIter.data(), do_dynamic_trim, dynamic_trim_tolerance, do_static_trim, static_trim_level);
			for ( ListIter<GraphArc *> arcIter(pathIter.data()) ; arcIter ; ++arcIter ) 
			{
				level++;
				a = arcIter.data();
				
				//cout << '(' << getKey(graphIdMap, a->source) << ' ' << getKey(graphIdMap, a->dest) << ' ' << a->weight << ')' << ' ';
				//cout << '(' << getKey(graphIdMap, a->dest) << ' ' << a->weight << ')' << endl;

				//printf("a->dest:%d\t", a->dest);
				//printf("getkey:%s\t", getKey(graphIdMap, a->dest).c_str());
				//printf("best_paths[%d] <- %s\n", node_pos, getKey(graphIdMap, a->dest).c_str());
				my_strcpy(best_paths[node_pos], getKey(graphIdMap, a->dest).c_str());
				node_pvalue[node_pos] = pvals[getKey(graphIdMap, a->dest)];
				node_pos++;
				//printf("(%s %lg)", getKey(graphIdMap, a->dest).c_str(), pvals[getKey(graphIdMap, a->dest)]);

				if (level > path_trim_level)
					continue;

				pathWeight += a->weight;
				addEdge(edge_list, getKey(graphIdMap, a->source), getKey(graphIdMap, a->dest), pvals[getKey(graphIdMap, a->dest)] - pvals[getKey(graphIdMap, a->source)]);
				min_edge = MIN(min_edge, pvals[getKey(graphIdMap, a->dest)] - pvals[getKey(graphIdMap, a->source)]);
				max_edge = MAX(max_edge, pvals[getKey(graphIdMap, a->dest)] - pvals[getKey(graphIdMap, a->source)]);
				//addNodeSize(node_size, getKey(graphIdMap, a->source), pval[getKey(graphIdMap, a->source)]);
				//addNodeSize(node_size, getKey(graphIdMap, a->dest), pval[getKey(graphIdMap, a->dest)] - pval[getKey(graphIdMap, a->source)]);
				node_size[getKey(graphIdMap, a->source)] = pvals[getKey(graphIdMap, a->source)];
				node_size[getKey(graphIdMap, a->dest)] = pvals[getKey(graphIdMap, a->dest)];
			}
			//printf("\n");
			//cout << pathWeight << '\n';
		}

		delete paths;
		delete[] graph.states;

		int node_it = 0;
		int node_attr_count = 3;
		for (map<string, double>::iterator it = node_size.begin(); it != node_size.end(); it++)
		{
			unsigned char scaled_value = (unsigned char)255 * (it->second - min_score) / (max_score - min_score);
			string color = "#";
			color += getHex(scaled_value) + "FF" + getHex(255 - scaled_value);
			//fprintf(stdout, "\"%s\" [color=\"%s\", style=filled];\n", signs[it->first].c_str(), color.c_str());
			my_strcpy(gnodes[node_it * node_attr_count + 0], it->first.c_str());
			my_strcpy(gnodes[node_it * node_attr_count + 1], signs[it->first].c_str());
			my_strcpy(gnodes[node_it * node_attr_count + 2], color.c_str());
			node_it++;
		}
		*node_count = node_it;

		int edge_it = 0;
		int edge_attr_count = 5;
		for (map<string, map<string, double> >::iterator it = edge_list.begin(); it != edge_list.end(); it++)
		{
			for (map<string, double>::iterator jt = it->second.begin(); jt != it->second.end(); jt++)
			{
				//float scaled_value = 0.75 + 15 * (jt->second - min_edge) / (max_edge - min_edge);
				//if (scaled_value < 0.5)
				//	scaled_value = 0.5;
				float scaled_value = jt->second;
				//fprintf(stdout, "%d\t \"%s\" -> \"%s\" [penwidth=%g, arrowhead=\"vee\", label=\"%s\", labelfontsize=14];\n", edge_it, signs[it->first].c_str(),
				//		signs[jt->first].c_str(), scaled_value, proteins[firstDiff(it->first, jt->first)]);
				char edge_name[strlen(it->first.c_str()) + strlen(jt->first.c_str()) + 2];
				sprintf(edge_name, "%s~%s", it->first.c_str(), jt->first.c_str());
				my_strcpy(gedges[edge_it * edge_attr_count + 0], edge_name);
				sprintf(gedges[edge_it * edge_attr_count + 1], "%g", scaled_value);

				string difference = diffPhenotype(it->first, jt->first);

				my_strcpy(gedges[edge_it * edge_attr_count + 2], difference.c_str());
				my_strcpy(gedges[edge_it * edge_attr_count + 3], "vee");
				my_strcpy(gedges[edge_it * edge_attr_count + 4], "14");
				edge_it++;
			}
		}
		*edge_count = edge_it;
	}
}

char **allocate(int x, int y)
{
	char **matrix = (char **)calloc(x, sizeof(char*));
	for(int i = 0; i < x; i++) 
	{
		matrix[i] = (char *)calloc(y, sizeof(char));
	}
	return matrix;
}

void destroy(char **v, int x)
{
	for (int i = 0; i < x; i++)
		free(v[i]);
	free(v);
}


/*
int main(int argc, char **argv)
{
	try
	{
		if (argc != 5)
		{
			fprintf(stderr, "usage: %s signs_file pvalue_file back_trace_start_string path_count\n", argv[0]);
			exit(-1);
		}
		char *pval_file_name = argv[2];
		char *signs_file_name = argv[1];
		char *back_trace_start = argv[3];
		int path_count = atoi(argv[4]);

		int trim_paths = 0;

		int data_count = 100000;
		int protein_count = strlen(back_trace_start);
		double *pvals = new double[data_count];
		int *s_values = new int[data_count * protein_count];
		char **signs = allocate(data_count, 255);
		int signs_row_num;
		char **proteins = allocate(10, 10);

		readPvals(signs_file_name, 
				pval_file_name, 
				pvals, 
				s_values, 
				signs, 
				signs_row_num,
				proteins);
	
		char **best_paths = allocate(path_count * protein_count, 255);;
		double node_pvalue[path_count * protein_count];
		char **gnodes = allocate(path_count * protein_count * 3, 255);
		char **gedges = allocate(path_count * protein_count * 5, 255);
		int node_count;
		int edge_count;
		c_analyze(signs, 
				s_values, 
				&signs_row_num, 
				&protein_count, 
				pvals, 
				proteins,
				&back_trace_start, 
				&path_count, 
				&trim_paths,
				//outputs
				best_paths, 
				node_pvalue,
				gnodes, 
				gedges,
				&node_count,
				&edge_count);
		//dumpGraphvizFile(graph, pvals, signs, proteins, graphIdMap, source, dest, path_count, "salam");
		//cout << graph;
		destroy(signs, data_count);
		destroy(proteins, 10);
		destroy(best_paths, path_count * protein_count);
		destroy(gnodes, path_count * protein_count * 3);
		destroy(gedges, path_count * protein_count * 5);
		delete[] pvals;
		delete[] s_values;
	}
	catch (Exception e)
	{
		fprintf(stderr, "%s\n", e.getMsg().c_str());
	}
	
	return 0;
}
*/
