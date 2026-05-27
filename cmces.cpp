#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <chrono>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <iostream>
#include <openbabel/obiter.h>
#include <utility>
#include <string>
#include <boost/graph/graphviz.hpp>
#include <fstream>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>
#include <boost/graph/connected_components.hpp>
#include <coroutine>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <fstream>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <deque>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <initializer_list>
#include <boost/graph/copy.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
 #include <string>
#include <algorithm>
#include <boost/graph/breadth_first_search.hpp>
#include <openbabel/bond.h>
#include <openbabel/math/vector3.h>
#include <openbabel/atom.h>
#include <unordered_map>


using namespace std;
using namespace boost;

using namespace OpenBabel;


typedef adjacency_list<vecS, vecS, undirectedS> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;

using Clique  = std::vector<Vertex>;
using Cliques = std::vector<Clique>;

std::vector<std::vector<Vertex>> all_cliques = std::vector<std::vector<Vertex>>();
std::vector<std::vector<Vertex>> all_maximal_subgraphs = std::vector<std::vector<Vertex>>(); 



std::vector<std::pair<Vertex, Vertex>> cartesian_product( Graph& g1, Graph& g2);
std::vector<struct ModularProductGraph> find_type1_connected_components(struct ModularProductGraph& graph);
struct MolecularGraph get_original_graph_from_modular_linegraph(struct ModularProductGraph& graph);
bool is_in(std::vector<Vertex> list, Vertex v);
std::vector<struct ModularProductGraph> alternate_find_type1_connected_components(struct ModularProductGraph& graph);
std::string to_mol(MolecularGraph &g);



int maximum = 0;
int uppr_bound = 0;
double time_for_cliques = 0;
int num_graphs = 0;
std::vector<int> clique_counter = std::vector<int>();
bool special_product = true; 
int  num_cliques_found = 0;
double path_finding_time = 0;
int edges_added = 0;
int edges_removed = 0; 




//Class for representing simple, labeled graphs
struct MolecularGraph{

    Graph g = Graph();

    std::map<Vertex, std::string> atoms =  std::map<Vertex, std::string>();
    std::map<std::pair<Vertex,Vertex>, std::string> bonds = std::map<std::pair<Vertex,Vertex>, std::string>();
    std::map<std::pair<Vertex,Vertex>, int> clique_map = std::map<std::pair<Vertex,Vertex>, int>();


    
  
    MolecularGraph() = default;



   
    MolecularGraph(const MolecularGraph& other)
        : g(other.g),                 
          atoms(other.atoms),         
          bonds(other.bonds),          
          clique_map(other.clique_map)
    {}


    MolecularGraph& operator=(const MolecularGraph& other) {
        if (this == &other) return *this;

       
        g = other.g;
        atoms = other.atoms;
        bonds = other.bonds;
        clique_map = other.clique_map;

        return *this;
    }
    


};
struct PathGraph{
    Graph g = Graph();
    std::map<Vertex, std::pair<Vertex,Vertex>> vertex_pair_map;
    std::map<std::pair<Vertex,Vertex>,Vertex> pair_vertex_map;



};


//Class for representing the line graph of MolecularGraph
struct LineGraph{

    Graph g = Graph();
    std::map<Vertex, std::string> atoms =  std::map<Vertex, std::string>();
    std::map<std::pair<Vertex,Vertex>, std::string> bonds = std::map<std::pair<Vertex,Vertex>, std::string>();
    MolecularGraph *org_graph = nullptr;
    std::map<Vertex, std::pair<Vertex, Vertex>> edge_map = std::map<Vertex, std::pair<Vertex, Vertex>>();
    std::map<Vertex, int> clique_map = std::map<Vertex, int>();
    int num = 0;
    std::vector<std::vector<std::vector<std::string>>> paths;
    std::vector<std::vector<int>> sizes;
    std::map<std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>, bool> map;
    bool dynamic = false;
    
    
    LineGraph() : org_graph(nullptr), dynamic(false) {}

    
    ~LineGraph() {
        if (dynamic) {
            delete org_graph;

        }
    }

    
    LineGraph(const LineGraph& other)
        : g(other.g), atoms(other.atoms), bonds(other.bonds),
          edge_map(other.edge_map), clique_map(other.clique_map),
          num(other.num), paths(other.paths), sizes(other.sizes),
          map(other.map), org_graph(nullptr), dynamic(other.dynamic) 
    {
        if (other.org_graph) {
            org_graph = new MolecularGraph(*other.org_graph);
            dynamic =   true;
        }
    }

    // Assignment operator (prevents memory leaks)
    LineGraph& operator=(const LineGraph& other) {
        if (this == &other) return *this; 

        
        if (dynamic) {
            delete org_graph;
        }

        // Copy data members
        g = other.g;
        atoms = other.atoms;
        bonds = other.bonds;
        edge_map = other.edge_map;
        clique_map = other.clique_map;
        num = other.num;
        paths = other.paths;
        sizes = other.sizes;
        map = other.map;

        // Deep copy for org_graph
        if (other.org_graph) {
            org_graph = new MolecularGraph(*other.org_graph);
            dynamic = true;  
        } else {
            org_graph = nullptr;
            dynamic = false;
        }

        return *this;
    }
};


//Class for representing the modular product of line graphs
struct ModularProductGraph{

    Graph g = Graph();
    std::map<Vertex, std::string> atoms =  std::map<Vertex, std::string>();
    std::map<std::pair<Vertex,Vertex>, std::string> bonds = std::map<std::pair<Vertex,Vertex>, std::string>();
    std::map<Vertex, std::pair<Vertex, Vertex>> node_map = std::map<Vertex, std::pair<Vertex, Vertex>>();
    struct LineGraph *g1;
    struct LineGraph *g2;
    std::map<Vertex, int> clique_map = std::map<Vertex, int>();



};





std::vector<std::vector<int>> all_pairs_shortest_paths(LineGraph& lg);
int max_size_type1 = 0;
struct MolecularGraph maximum_graph = MolecularGraph();
std::vector<MolecularGraph> all_maximum_graphs = std::vector<MolecularGraph>();
std::vector<std::vector<std::vector<Vertex>>> clique_vectors = std::vector<std::vector<std::vector<Vertex>>>();
std::map<LineGraph, int> clique_vector_map;
std::vector<std::map<std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>, bool>> all_maps = std::vector<std::map<std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>, bool>>();
std::vector<PathGraph> pgraphs = std::vector<PathGraph>();
std::vector<Vertex> refine_graph(ModularProductGraph &graph);

bool equivalent_paths(Vertex u1, Vertex u2, LineGraph &g1, Vertex v1, Vertex v2, LineGraph &g2);

struct Collector {
    Cliques& target;

    void clique(auto const& clique, Graph const&) const {
        for (auto& t = target.emplace_back(); Graph::vertex_descriptor v : clique)
            t.push_back(v);
    }
};








//Class used for "yielding" cliques in the clique finding function, similar to the keyword yield in python
class CliqueGenerator {
public:
    struct promise_type {
        std::set<Vertex> current_clique;

        auto yield_value(std::set<Vertex> clique) {
            current_clique = std::move(clique);
            return std::suspend_always{};
        }

        auto return_void() { return std::suspend_always{}; }

        static auto get_return_object_on_allocation_failure() { return CliqueGenerator{nullptr}; }
        auto get_return_object() { return CliqueGenerator{std::coroutine_handle<promise_type>::from_promise(*this)}; }

        auto initial_suspend() { return std::suspend_always{}; }
        auto final_suspend() noexcept { return std::suspend_always{}; }

        void unhandled_exception() {
            std::terminate(); 
        }
    };

    using handle_type = std::coroutine_handle<promise_type>;

    CliqueGenerator(handle_type h) : handle_(h) {}

    ~CliqueGenerator() {
        if (handle_) handle_.destroy();
    }

    bool move_next() {
        handle_.resume();
        return !handle_.done();
    }

    std::set<Vertex> current_value() {
        return handle_.promise().current_clique;
    }

private:
    handle_type handle_;
};


CliqueGenerator bron_kerbosch_driver(ModularProductGraph &graph);



//Class for yielding maximal common subgraphs, similar to the yield keyword in python
class MCSGenerator {
public:
    struct promise_type {
        MolecularGraph current_subgraph;

        auto yield_value(MolecularGraph subgraph) {
            current_subgraph = std::move(subgraph);
            return std::suspend_always{};
        }

        auto return_void() { return std::suspend_always{}; }

        static auto get_return_object_on_allocation_failure() { return MCSGenerator{nullptr}; }
        auto get_return_object() { return MCSGenerator{std::coroutine_handle<promise_type>::from_promise(*this)}; }

        auto initial_suspend() { return std::suspend_always{}; }
        auto final_suspend() noexcept { return std::suspend_always{}; }

        void unhandled_exception() {
            std::terminate(); 
        }
    };

    using handle_type = std::coroutine_handle<promise_type>;

    MCSGenerator(handle_type h) : handle_(h) {}

    ~MCSGenerator() {
        if (handle_) handle_.destroy();
    }

    bool move_next() {
        handle_.resume();
        return !handle_.done();
    }

    MolecularGraph current_value() {
        return handle_.promise().current_subgraph;
    }

private:
    handle_type handle_;
};








//Constructs the modular product of two line graphs, WITHOUT type-0 edges, i.e. this function constructs the direct product of two line graphs
struct ModularProductGraph modular_graph_product(struct LineGraph& g1, struct LineGraph& g2){

    struct ModularProductGraph modular_product;


    modular_product.g1 = &g1;
    modular_product.g2 = &g2;



    for(auto &node: cartesian_product(g1.g, g2.g)){



        if(g1.atoms[node.first]!=g2.atoms[node.second]){
            continue;
        }

        
        auto edge1 = g1.edge_map[node.first];
        auto edge2 = g2.edge_map[node.second];

       


        if( !( (g1.org_graph->atoms[edge1.first]==g2.org_graph->atoms[edge2.first] && g1.org_graph->atoms[edge1.second]==g2.org_graph->atoms[edge2.second]) ||  (g1.org_graph->atoms[edge1.first]==g2.org_graph->atoms[edge2.second] && g1.org_graph->atoms[edge1.second]==g2.org_graph->atoms[edge2.first]))){
            continue;
        }
         

        auto v = add_vertex(modular_product.g);
        modular_product.node_map[v] = node;
        modular_product.atoms[v] = g1.atoms[node.first];
        modular_product.clique_map[v] = g1.clique_map[node.first];

    }
    
    for (const auto& node1 : boost::make_iterator_range(boost::vertices(modular_product.g))){

        for (const auto& node2 : boost::make_iterator_range(boost::vertices(modular_product.g))){
            if(node1<=node2){
                continue;
            }

            auto tup1 = modular_product.node_map[node1];
            auto tup2 = modular_product.node_map[node2];

            auto edge1 = g1.edge_map[tup1.first];
            auto edge2 = g1.edge_map[tup2.first];

        

            if(tup1.first==tup2.first || tup1.second==tup2.second){
                continue;
            }
            else if(false && !edge(tup1.first, tup2.first, g1.g).second and !edge(tup1.second, tup2.second, g2.g).second ){
                auto e = add_edge(node1, node2, modular_product.g);
                auto pair1 = std::pair<Vertex, Vertex>(node1, node2);
                auto pair2 = std::pair<Vertex, Vertex>(node2, node1);
                modular_product.bonds[pair1] = "NONE";
                modular_product.bonds[pair2] = "NONE";

            }
            else if(edge(tup1.first, tup2.first, g1.g).second && edge(tup1.second, tup2.second, g2.g).second && g1.bonds[std::pair<Vertex, Vertex>(tup1.first, tup2.first)] == g2.bonds[std::pair<Vertex, Vertex>(tup1.second, tup2.second)]){
                auto e = add_edge(node1, node2, modular_product.g);
                auto pair1 = std::pair<Vertex, Vertex>(node1, node2);
                auto pair2 = std::pair<Vertex, Vertex>(node2, node1);
                modular_product.bonds[pair1] = g1.bonds[std::pair<Vertex, Vertex>(tup1.first, tup2.first)];
                modular_product.bonds[pair2] = g1.bonds[std::pair<Vertex, Vertex>(tup1.first, tup2.first)];



            }

        }

    }


    return modular_product;



}



//Takes a simple labeled graph(representing a molecule) and returns its line graph
LineGraph molecular_graph_to_line_graph(MolecularGraph &graph){

    LineGraph line_graph;

    auto pointer = &graph;
    line_graph.org_graph = pointer;


    

    

    for (const auto& edge : boost::make_iterator_range(boost::edges(graph.g))) {
        auto v = add_vertex(line_graph.g);
        auto edge_pair = std::pair<Vertex, Vertex>(edge.m_target, edge.m_source);
        line_graph.atoms[v] = graph.bonds[edge_pair];
        line_graph.edge_map[v] = edge_pair;
        line_graph.clique_map[v] = graph.clique_map[edge_pair];
    }

     for (const auto& node1 : boost::make_iterator_range(boost::vertices(line_graph.g))){

        for (const auto& node2 : boost::make_iterator_range(boost::vertices(line_graph.g))){
            if(node1<=node2){
                continue;
            }
            


            auto edge1 = line_graph.edge_map[node1];
            auto edge2 = line_graph.edge_map[node2];

            auto pair1 = std::pair<Vertex, Vertex>(node1, node2);
            auto pair2 = std::pair<Vertex, Vertex>(node2, node1);

            if(edge1.first==edge2.first){
                auto e = add_edge(node1, node2, line_graph.g);
                line_graph.bonds[pair1] = graph.atoms[edge1.first];
                line_graph.bonds[pair2] = graph.atoms[edge1.first];
            }
            else if(edge1.first==edge2.second){
                auto e = add_edge(node1, node2, line_graph.g);
                line_graph.bonds[pair1] = graph.atoms[edge1.first];
                line_graph.bonds[pair2] = graph.atoms[edge1.first];

            }
            else if(edge1.second==edge2.second){
                auto e = add_edge(node1, node2, line_graph.g);
                line_graph.bonds[pair1] = graph.atoms[edge1.second];
                line_graph.bonds[pair2] = graph.atoms[edge1.second];


            }
            else if(edge1.second==edge2.first){
                auto e = add_edge(node1, node2, line_graph.g);
                line_graph.bonds[pair1] = graph.atoms[edge1.second];
                line_graph.bonds[pair2] = graph.atoms[edge1.second];

            }

        }

     }


    std::vector<std::vector<std::vector<std::string>>> paths;
   
    line_graph.paths = paths;



    return line_graph;


}


//helper functions for checking inclusion in various ways:
/*****************************************/
bool is_in(std::set<std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>> list, std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>> v){
    for(const auto& u: list){
        if(u==v){
            return true;
        }
    }
    return false;
}

bool is_in(std::set<std::pair<Vertex, Vertex>> list, std::pair<Vertex, Vertex> v){
    for(const auto& u: list){
        if(u==v){
            return true;
        }
    }
    return false;
}

bool is_in(std::vector<Vertex> list, Vertex v){
    for(const auto& u: list){
        if(u==v){
            return true;
        }
    }
    return false;
}

bool is_in(std::vector<int> list, int v){
    for(const auto& u: list){
        if(u==v){
            return true;
        }
    }
    return false;
}

bool is_in(std::set<Vertex> list, Vertex v){
    for(const auto& u: list){
        if(u==v){
            return true;
        }
    }
    return false;
}

bool is_in(std::vector<std::pair<Vertex, Vertex>> list, std::pair<Vertex, Vertex> v){
    for(auto& u: list){

        if(u.first==v.first && u.second == v.second){
            return true;
        }
    }
    return false;
}

bool is_in2(std::vector<std::pair<Vertex, Vertex>> list, std::pair<Vertex, Vertex> v){
    for(auto& u: list){

        if(u.first==v.first && u.second == v.second){
            return true;
        }
    }
    return false;
}

bool is_in(std::set<Vertex> triangle, std::vector<std::set<Vertex>> seen_triangles){

    for(auto tri : seen_triangles){
        if(tri==triangle){
            return true;
        }
    }

    return false;
}
//helper functions for checking inclusion in various ways:
/**************************************************/




//Returns the cartesian product of an arbitrary number of sets
std::vector<std::set<Vertex>> cartesian(std::vector<std::vector<Vertex>>& sets) {
    std::vector<std::set<Vertex>> temp(1, std::set<Vertex>());
    for (int i = 0; i < sets.size(); i++) {
        std::vector<std::set<Vertex>> newTemp;
        for (const std::set<Vertex>& product : temp) {
            for (const Vertex& element : sets[i]) {
                std::set<Vertex> tempCopy = product;
                tempCopy.insert(element);
                newTemp.push_back(tempCopy);
            }
        }
        temp = newTemp;
    }

    return temp;

}




//Takes a modular product graph and a list of vertices, and deletes the vertices from the graph
ModularProductGraph delete_set_of_vertices(ModularProductGraph &graph, std::set<Vertex> &delete_vertices){


    std::set<Vertex> good_vertices;

    for (const auto& vertex : boost::make_iterator_range(boost::vertices(graph.g))){
        if(is_in(delete_vertices, vertex)){
            continue;
        }
        good_vertices.insert(vertex);
    }

    struct ModularProductGraph subgraph;

    auto vertex_map = std::map<Vertex, Vertex>();

    subgraph.g1 = graph.g1;
    subgraph.g2 = graph.g2;

    for(const auto& vertex: good_vertices){

        auto v = add_vertex(subgraph.g);
        vertex_map[vertex] = v;
        subgraph.atoms[v] = graph.atoms[vertex];
        subgraph.node_map[v] = graph.node_map[vertex];
        subgraph.clique_map[v] = graph.g1->clique_map[graph.node_map[vertex].first]; 

    }

    for (const auto& edge : boost::make_iterator_range(boost::edges(graph.g))){
        if(good_vertices.find(edge.m_source)!=good_vertices.end() && good_vertices.find(edge.m_target)!=good_vertices.end() && graph.bonds[std::pair<Vertex, Vertex>(edge.m_target, edge.m_source)] != "NONE"){//(is_in(clique, edge.m_source) && is_in(clique, edge.m_target) && graph.bonds[std::pair<Vertex, Vertex>(edge.m_target, edge.m_source)] != "NONE"){
            auto e = add_edge(vertex_map[edge.m_source], vertex_map[edge.m_target], subgraph.g);
            subgraph.bonds[std::pair<Vertex, Vertex>(vertex_map[edge.m_target],vertex_map[edge.m_source])] = graph.bonds[std::pair<Vertex, Vertex>(edge.m_target, edge.m_source)];
            subgraph.bonds[std::pair<Vertex, Vertex>(vertex_map[edge.m_source],vertex_map[edge.m_target])] = graph.bonds[std::pair<Vertex, Vertex>(edge.m_target, edge.m_source)];
        }

    }

    return subgraph;






}


//Returns an induced subgraph of a modular product graph containing the vertices in clique. Also checks if there is a bad triangle
//in the given clique, in which case it returns three representative modular product graphs, as outlined in the paper
std::vector<ModularProductGraph> induced_subgraph(const std::set<Vertex> &clique, struct ModularProductGraph& graph){

    struct ModularProductGraph subgraph;

    auto vertex_map = std::map<Vertex, Vertex>();

    subgraph.g1 = graph.g1;
    subgraph.g2 = graph.g2;

    for(const auto& vertex: clique){

        auto v = add_vertex(subgraph.g);
        vertex_map[vertex] = v;
        subgraph.atoms[v] = graph.atoms[vertex];
        subgraph.node_map[v] = graph.node_map[vertex];
        subgraph.clique_map[v] = graph.g1->clique_map[graph.node_map[vertex].first]; 

    }

    for (const auto& edge : boost::make_iterator_range(boost::edges(graph.g))){
        if(clique.find(edge.m_source)!=clique.end() && clique.find(edge.m_target)!=clique.end() && graph.bonds[std::pair<Vertex, Vertex>(edge.m_target, edge.m_source)] != "NONE"){//(is_in(clique, edge.m_source) && is_in(clique, edge.m_target) && graph.bonds[std::pair<Vertex, Vertex>(edge.m_target, edge.m_source)] != "NONE"){
            auto e = add_edge(vertex_map[edge.m_source], vertex_map[edge.m_target], subgraph.g);
            subgraph.bonds[std::pair<Vertex, Vertex>(vertex_map[edge.m_target],vertex_map[edge.m_source])] = graph.bonds[std::pair<Vertex, Vertex>(edge.m_target, edge.m_source)];
            subgraph.bonds[std::pair<Vertex, Vertex>(vertex_map[edge.m_source],vertex_map[edge.m_target])] = graph.bonds[std::pair<Vertex, Vertex>(edge.m_target, edge.m_source)];
        }

    }

    std::vector<std::vector<Vertex>> bad_triangles;

    std::vector<std::set<Vertex>> seen_triangles;

    if(clique.size()!=3){
        return {subgraph};
    }
    
 

    for (const auto& e : boost::make_iterator_range(boost::edges(subgraph.g))){
        if(subgraph.bonds[std::pair<Vertex, Vertex>(vertex_map[e.m_target],vertex_map[e.m_source])]=="NONE"){
            continue;
        }
        for (const auto& vertex : boost::make_iterator_range(boost::vertices(subgraph.g))){
            if(edge(vertex, e.m_target, subgraph.g).second and edge(vertex, e.m_source, subgraph.g).second and subgraph.bonds[std::pair<Vertex, Vertex>(vertex,e.m_source)]!="NONE" and subgraph.bonds[std::pair<Vertex, Vertex>(vertex,e.m_target)]!="NONE"){
                std::set<Vertex> triangle1;
                std::set<Vertex> triangle2;

                triangle1.insert(subgraph.g1->edge_map[subgraph.node_map[vertex].first].first);
                triangle1.insert(subgraph.g1->edge_map[subgraph.node_map[vertex].first].second);

                triangle1.insert(subgraph.g1->edge_map[subgraph.node_map[e.m_source].first].first);
                triangle1.insert(subgraph.g1->edge_map[subgraph.node_map[e.m_source].first].second);

                triangle1.insert(subgraph.g1->edge_map[subgraph.node_map[e.m_target].first].first);
                triangle1.insert(subgraph.g1->edge_map[subgraph.node_map[e.m_target].first].second);


                triangle2.insert(subgraph.g2->edge_map[subgraph.node_map[vertex].second].first);
                triangle2.insert(subgraph.g2->edge_map[subgraph.node_map[vertex].second].second);

                triangle2.insert(subgraph.g2->edge_map[subgraph.node_map[e.m_source].second].first);
                triangle2.insert(subgraph.g2->edge_map[subgraph.node_map[e.m_source].second].second);

                triangle2.insert(subgraph.g2->edge_map[subgraph.node_map[e.m_target].second].first);
                triangle2.insert(subgraph.g2->edge_map[subgraph.node_map[e.m_target].second].second); 

                if(triangle1.size()!=triangle2.size()){
                    std::vector<Vertex> bad_triangle = {vertex, e.m_target, e.m_source};
                    std::set<Vertex> triangle = {vertex, e.m_target, e.m_source};
                    if(is_in(triangle, seen_triangles)){
                        continue;
                    }
                    bad_triangles.push_back(bad_triangle);
                    seen_triangles.push_back(triangle);
                }
            }
        }

    }

    if(bad_triangles.size()==0){
        return {subgraph};
    }

    std::vector<ModularProductGraph> results;


    
    auto remove_vertex_sets = cartesian(bad_triangles);

    std::vector<std::set<Vertex>> final_remove_vertex_sets;

    std::vector<bool> keep_set;

    for(int i=0; i<remove_vertex_sets.size(); i++){
        keep_set.push_back(true);
    }

    for(int i = 0; i<remove_vertex_sets.size(); i++){
        for(int j=0; j<remove_vertex_sets.size(); j++){
            if(i==j){
                continue;
            }
            if(std::includes(remove_vertex_sets[i].begin(), remove_vertex_sets[i].end(), remove_vertex_sets[j].begin(), remove_vertex_sets[j].end() )){
                keep_set[j]=false;
            }


        }

    }

    for(int i=0; i<remove_vertex_sets.size(); i++){
        if(keep_set[i]){
            final_remove_vertex_sets.push_back(remove_vertex_sets[i]);
        }
    }

    for(auto remove_set : final_remove_vertex_sets){

        auto curr_subgraph = delete_set_of_vertices(subgraph, remove_set);

        results.push_back(curr_subgraph);


    }





    return results;
    



}



//Constructs the original simple, labeled graph given a line graph
MolecularGraph get_original_graph_from_linegraph(LineGraph& g1, LineGraph& g2, std::vector<std::pair<Vertex,Vertex>> node_pairs){


    struct MolecularGraph orig_graph;


    std::set<Vertex> vertices;

    std::set<std::pair<Vertex, Vertex>> edges;


    for (auto node_pair : node_pairs){
        auto v = node_pair.first;

        if(vertices.find(g1.edge_map[v].first)==vertices.end()){

            vertices.insert(g1.edge_map[v].first);
        }

        if(vertices.find(g1.edge_map[v].second)==vertices.end()){
            vertices.insert(g1.edge_map[v].second);
        }
        

        edges.insert(g1.edge_map[v]);
    }



    std::map<Vertex, Vertex> vertex_map;

    for(auto vertex : vertices){
        auto v = add_vertex(orig_graph.g);
        vertex_map[vertex] = v;
        orig_graph.atoms[v] = g1.org_graph->atoms[vertex];


    }
    for (auto edge : boost::make_iterator_range(boost::edges(g1.org_graph->g))){
        
    
        if(edges.find(std::pair<Vertex, Vertex>(edge.m_source, edge.m_target))!=edges.end()  or edges.find(std::pair<Vertex, Vertex>(edge.m_target, edge.m_source))!=edges.end()){
           
            auto e = add_edge(vertex_map[edge.m_source], vertex_map[edge.m_target], orig_graph.g);
            orig_graph.bonds[std::pair<Vertex, Vertex>(vertex_map[edge.m_source], vertex_map[edge.m_target])] = g1.org_graph->bonds[(std::pair<Vertex, Vertex>(edge.m_source, edge.m_target))];
            orig_graph.bonds[std::pair<Vertex, Vertex>(vertex_map[edge.m_target], vertex_map[edge.m_source])] = g1.org_graph->bonds[(std::pair<Vertex, Vertex>(edge.m_source, edge.m_target))];

            orig_graph.clique_map[std::pair<Vertex, Vertex>(vertex_map[edge.m_target], vertex_map[edge.m_source])] = g1.org_graph->clique_map[(std::pair<Vertex, Vertex>(edge.m_source, edge.m_target))];
            orig_graph.clique_map[std::pair<Vertex, Vertex>(vertex_map[edge.m_source], vertex_map[edge.m_target])] = g1.org_graph->clique_map[(std::pair<Vertex, Vertex>(edge.m_source, edge.m_target))];
           
        }



    }


    return orig_graph;



}

//possible implementation of alternative maximal common subgraph algorithm. DOES NOT WORK YET
MCSGenerator alt_maximal_common_subgraphs(struct LineGraph &g1, struct LineGraph &g2){

    std::vector<std::vector<std::pair<Vertex,Vertex>>> stack;
    
    
     for (const auto& v1 : boost::make_iterator_range(boost::vertices(g1.g))){
    
        for (const auto& v2 : boost::make_iterator_range(boost::vertices(g2.g))){
            if(g1.atoms[v1]==g2.atoms[v2]){
                std::pair<Vertex,Vertex> first_pair;
                first_pair.first = v1;
                first_pair.second = v2;
                std::vector<std::pair<Vertex,Vertex>> first = {first_pair};


                stack.push_back(first);
            }

        }

     }

     int max = 0;
     

     std::vector<std::pair<Vertex,Vertex>> bad_pairs;

     while(stack.size()>0){

        auto curr = stack.back();

        stack.pop_back();

        if(curr.size()>max){
            max = curr.size();
        }

        bool is_maximal = true;

        for(auto pair : curr){
            bool bad = false;
            for(auto p : bad_pairs){
                if(pair.first==p.first and pair.second == p.second){
                    bad = true;
                    break;
                }
            }
            if(bad){
                continue;
            }
            bool good = false;
            auto v1 = pair.first;
            auto v2 = pair.second;
            for(const auto new_v1 : make_iterator_range(adjacent_vertices(v1, g1.g))){
                for(const auto new_v2 : make_iterator_range(adjacent_vertices(v2, g2.g))){
                    
                    if(g1.atoms[new_v1]==g2.atoms[new_v2] and g1.bonds[std::pair(v1,new_v1)]==g2.bonds[std::pair(v2, new_v2)]){
                        


                        bool is_present = false;
                        
                        good = true;
                        for(auto p : curr){
                            if((p.first==new_v1 or p.second==new_v2)){
                                is_present = true;
                                break;
                            }
                            if(g1.bonds[std::pair(p.first, new_v1)]!=g2.bonds[std::pair(p.second, new_v2)]){
                                is_present = true;
                                break;
                            }
                        }
                        
                        if(is_present){
                            continue;
                        }
                        
                        is_maximal = false;
                        auto new_item = curr;
                        new_item.push_back(std::pair(new_v1, new_v2));
                        stack.push_back(new_item);
                        

                    }
                }
            }
            if(not good){
                bad_pairs.push_back(pair);
            }
        }
        if(is_maximal){

            auto orig = get_original_graph_from_linegraph(g1, g2, curr);
             

            co_yield orig;

        }

     }

     



     co_return;



}




//Returns a generator which yields all maximal common subgraphs between g1 and g2
//If global var special_product is true it removes certain type-0 edges from the modular product of g1 and g2
//and otherwise it doesn't 
MCSGenerator maximal_common_subgraphs(struct LineGraph &g1, struct LineGraph &g2){

    std::vector<Vertex> clique_set_prime;
    std::vector<Vertex> clique_set2_prime;

    for(auto vertex :  boost::make_iterator_range(boost::vertices(g1.g))){
        clique_set_prime.push_back(g1.clique_map[vertex]);
        
    }

    std::sort(clique_set_prime.begin(), clique_set_prime.end());
    std::sort(clique_set2_prime.begin(), clique_set2_prime.end());

    bool bad_clique = false;

    for(auto &test_clique : clique_vectors[0]){
        if(std::includes(test_clique.begin(), test_clique.end(), clique_set_prime.begin(), clique_set_prime.end())){
            bad_clique = true;
            break;
        }
        

    }
    if(bad_clique){
        co_return;
    }
    
    

    int max = 0;
    struct ModularProductGraph max_graph;

    

    auto mod_product = modular_graph_product(g1, g2);

    std::vector<ModularProductGraph> subgraphs;

    if(special_product){
        subgraphs =alternate_find_type1_connected_components(mod_product);
    }
    else{
        subgraphs =find_type1_connected_components(mod_product);
    }

    
  
    std::vector<MolecularGraph> results;

    for( auto& subgraph: subgraphs){


        if(num_vertices(subgraph.g)<maximum){
            continue;
        }

        if(num_vertices(g1.g)<maximum){
            co_return;
        }

        
        
        auto cliques =  bron_kerbosch_driver(subgraph);
        while (cliques.move_next()){

            auto clique = cliques.current_value();

            num_cliques_found++;
            

            if(clique.size()<maximum){
                co_return;
            }

            

            std::vector<Vertex> clique_set;
            std::vector<Vertex> clique_set2;

            for(auto &vertex : clique){
                clique_set.push_back(subgraph.clique_map[vertex]);
                clique_set2.push_back(subgraph.g2->clique_map[subgraph.node_map[vertex].second]);
            }

            std::sort(clique_set.begin(), clique_set.end());
            std::sort(clique_set2.begin(), clique_set2.end());

            bool bad_clique = false;

            for(auto &test_clique : clique_vectors[0]){
                if(std::includes(test_clique.begin(), test_clique.end(), clique_set.begin(), clique_set.end())){
                    bad_clique = true;
                    break;
                }
                

            }
            if(bad_clique){
                clique_vectors[g2.num].push_back(clique_set2);
                continue;
            }
            for(auto &test_clique : clique_vectors[g2.num]){
                if(std::includes(test_clique.begin(), test_clique.end(), clique_set2.begin(), clique_set2.end())){ //or std::includes(clique_set2.begin(), clique_set2.end(), test_clique.begin(), test_clique.end())){
                    bad_clique = true;
                    break;
                }
                

            }

            if(bad_clique){
                clique_vectors[0].push_back(clique_set);
                continue;
            }
            



            



            

            auto ind_subgraph = induced_subgraph(clique, subgraph);

            //std::vector<ModularProductGraph> smaller_subgraphs;
            //smaller_subgraphs = {ind_subgraph};

            
            
            for(auto& small_subgraph: ind_subgraph){

                

                if(num_vertices(small_subgraph.g)>=maximum){






                    auto molecular_graph = get_original_graph_from_modular_linegraph(small_subgraph);

                    results.push_back(molecular_graph);

                    clique_counter[g2.num-1]++;

                    co_yield molecular_graph;

                                      

                }
                
            }
            clique_vectors[0].push_back(clique_set);
            clique_vectors[g2.num].push_back(clique_set2);

            int size = 0;

            for(int i=0; i<clique_vectors.size(); i++){
                size = size + clique_vectors[i].size();
            }
            if( size>100000000){
                std::vector<std::vector<std::vector<Vertex>>> outer_vector;

                for(int i=0; i<num_graphs; i++){
                    std::vector<std::vector<Vertex>> inner_vector;
                    outer_vector.push_back(inner_vector);
                    
                }
                clique_vectors = outer_vector;
            }
            
        }

            
            
        }


    }








//Returns the cartesian product of two sets
std::vector<std::pair<Vertex, Vertex>> cartesian_product( Graph& g1, Graph& g2){

    std::vector<std::pair<Vertex, Vertex>> product;

    for (const auto& vi : boost::make_iterator_range(boost::vertices(g1))){
        for (const auto& vu : boost::make_iterator_range(boost::vertices(g2))){

            auto pair = std::pair<Vertex, Vertex>(vi,vu);
            product.push_back(pair);
        }
    }

    return product;

}




struct StackHelper{
    Vertex u_prime;
    Vertex v_prime;
    std::set<Vertex> seen_u_vertices;
    std::set<Vertex> seen_v_vertices;
    std::vector<std::pair<Vertex, Vertex>> path; 
};



//Performs the isomorphic path finding heuristic outlined in the paper. 
//Performs an exhaustive simultaneous breadth first search in g1 and g2 starting in u1 and v1
//extending isomorphic paths in each step. Returns all pairs reachable from u1 and v1.
//If number of steps exceeds 10000 it terminates and simply returns all pairs, thereby pruning no type-0 edges.
std::vector<std::pair<Vertex, Vertex>> reachable_pairs(Vertex u1, LineGraph &g1, Vertex v1, LineGraph &g2, std::vector<Vertex> good_vertices1, std::vector<Vertex> good_vertices2, int max_num_edges){


    std::vector<std::pair<Vertex, Vertex>> all_pairs;
    for (auto ud : make_iterator_range(vertices(g1.g))){
           
        for (auto vd : make_iterator_range(vertices(g2.g))){
            all_pairs.push_back(std::pair(ud,vd));

        }
    }


    
    std::vector<std::pair<Vertex, Vertex>> reachable;
    std::deque<StackHelper> stack;

    StackHelper first;
    first.u_prime = u1;
    first.v_prime = v1;

    std::set<std::pair<Vertex,Vertex>> seen_pairs;
    
    first.seen_u_vertices.insert(u1);
    first.seen_v_vertices.insert(v1);
    stack.push_back(first);

    int iter = 0;

    int pop_iter = 0;

    auto quad = first;

    while(stack.size()>0){
        iter++;
        

        if(iter>10000){
            return all_pairs;
        }

    
        quad = stack.back();
        
        auto u_prime = quad.u_prime;
        auto v_prime = quad.v_prime;

        stack.pop_back();
        
        bool abort = false;

        for(auto stack_item : stack){
            if(stack_item.u_prime==u_prime and stack_item.v_prime==v_prime){
                if(std::includes(quad.seen_u_vertices.begin(), quad.seen_u_vertices.end(), stack_item.seen_u_vertices.begin(), stack_item.seen_u_vertices.end()) and std::includes(quad.seen_v_vertices.begin(), quad.seen_v_vertices.end(), stack_item.seen_v_vertices.begin(), stack_item.seen_v_vertices.end()) ){
                    abort = true;
                    break;
                }
            }
        }

        if(abort){
            continue;
        }
        
        for (auto ud : make_iterator_range(adjacent_vertices(u_prime, g1.g))){
           
            if(quad.seen_u_vertices.find(ud)!=quad.seen_u_vertices.end()){
                continue;
            }
            if(!std::binary_search(good_vertices1.begin(), good_vertices1.end(), ud)){
                    continue;
                }
            for (auto vd : make_iterator_range(adjacent_vertices(v_prime, g2.g))){

                auto edge1 = g1.edge_map[ud];
                auto edge2 = g2.edge_map[vd];

                if(!( (g1.org_graph->atoms[edge1.first]==g2.org_graph->atoms[edge2.first] && g1.org_graph->atoms[edge1.second]==g2.org_graph->atoms[edge2.second]) ||  (g1.org_graph->atoms[edge1.first]==g2.org_graph->atoms[edge2.second] && g1.org_graph->atoms[edge1.second]==g2.org_graph->atoms[edge2.first]))){
                    continue;
                }

            
                if(g1.atoms[ud]!=g2.atoms[vd]){
                    continue;
                }
                if(g1.bonds[std::pair<Vertex, Vertex>(u_prime, ud)]!=g2.bonds[std::pair<Vertex, Vertex>(v_prime, vd)]){
                    continue;
                }
        
                if(!std::binary_search(good_vertices2.begin(), good_vertices2.end(), vd)){
                    continue;
                }
                if(quad.seen_v_vertices.find(vd)!=quad.seen_v_vertices.end()){
                    continue;
                }
                else{

                    if(!is_in(reachable, std::pair(ud,vd))){
                        reachable.push_back(std::pair(ud,vd));
                    }
                    
                    StackHelper new_quad;
                    new_quad.u_prime = ud;
                    new_quad.v_prime = vd;

                    for(auto h : quad.seen_u_vertices){
                        new_quad.seen_u_vertices.insert(h);
                    }
                    for(auto h : quad.seen_v_vertices){
                        new_quad.seen_v_vertices.insert(h);
                    }

                    new_quad.seen_u_vertices.insert(ud);
                    new_quad.seen_v_vertices.insert(vd);

                    stack.push_back(new_quad);

                }
            }


        }
        
        
    } 

    return reachable;
    

    
}











//Finds the type-A connected components of a modular product graph, without the removal of certain type-0 edges
std::vector<struct ModularProductGraph> find_type1_connected_components(struct ModularProductGraph& graph){

    std::vector<std::pair<Vertex, Vertex>> vertex_pairs;

    
    
    std::vector<struct ModularProductGraph> results;


    std::vector<int> component (boost::num_vertices (graph.g));
    int num_components = boost::connected_components (graph.g, &component[0]);



    

    std::map<int, std::vector<Vertex>> component_map = std::map<int, std::vector<Vertex>>();


    for(int i=0; i<num_components; i++){
        component_map[i] = std::vector<Vertex>();
    }

    for (const auto& vi : boost::make_iterator_range(boost::vertices(graph.g))){
        component_map[component[vi]].push_back(vi);

    }

    auto g1 = graph.g1;
    auto g2 = graph.g2;
    

    std::map<Vertex, Vertex> inverse_clique_map;

    for (const auto& vi : boost::make_iterator_range(boost::vertices(g1->g))){
        inverse_clique_map[g1->clique_map[vi]] = vi;
    }



    for(int i=0; i<num_components; i++){



        auto curr_vertices = component_map[i];

        if(curr_vertices.size()==1 or curr_vertices.size()<maximum){
            continue;
        }

        struct ModularProductGraph new_mod_product;

        new_mod_product.g1 = graph.g1;
        new_mod_product.g2 = graph.g2;

        std::vector<Vertex> good_vertices1;
        std::vector<Vertex> good_vertices2;

        std::map<std::pair<Vertex, Vertex>, Vertex> pair_vertex_map;

        for(auto vertex: curr_vertices){
            
            auto v = add_vertex(new_mod_product.g);
            new_mod_product.node_map[v] = graph.node_map[vertex];
            new_mod_product.atoms[v] = graph.atoms[vertex];
            new_mod_product.clique_map[v] = graph.clique_map[vertex];
            pair_vertex_map[graph.node_map[vertex]] = v;
            good_vertices1.push_back(graph.node_map[vertex].first);
            good_vertices2.push_back(graph.node_map[vertex].second);
            
        }

        std::sort(good_vertices1.begin(), good_vertices1.end());
        std::sort(good_vertices2.begin(), good_vertices2.end());

        

        double num_iter = (num_vertices(new_mod_product.g) * num_vertices(new_mod_product.g))/2;

        int curr_num = 0;
        

        std::set<std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>> pairs;

        std::map<std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>, bool> pair_map;

        int iter = 0;
        
        for (const auto& node1 : boost::make_iterator_range(boost::vertices(new_mod_product.g))){
            for (const auto& node2 : boost::make_iterator_range(boost::vertices(new_mod_product.g))){



                if(node1<=node2){
                    continue;
                }

                iter++;
                
                
                if(edge(node1, node2, new_mod_product.g).second){
                    continue;
                }
                curr_num ++;
                
                

                auto tup1 = new_mod_product.node_map[node1];
                auto tup2 = new_mod_product.node_map[node2];

                auto edge1 = g1->edge_map[tup1.first];
                auto edge2 = g1->edge_map[tup2.first];




                
               
                
                
                if(tup1.first==tup2.first || tup1.second==tup2.second){
                    continue;
                }



                

                auto pair_pair1 = std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>(std::pair<Vertex, Vertex>(graph.g1->clique_map[tup1.first], tup1.second), std::pair<Vertex, Vertex>(graph.g1->clique_map[tup2.first], tup2.second));
                auto pair_pair2 = std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>(std::pair<Vertex, Vertex>(graph.g1->clique_map[tup2.first], tup2.second), std::pair<Vertex, Vertex>(graph.g1->clique_map[tup1.first], tup1.second));

                
                bool condition_two = false;

                std::vector<std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex,Vertex>>> new_edges;

                if(!edge(tup1.first, tup2.first, g1->g).second && !edge(tup1.second, tup2.second, g2->g).second and (pair_map[(std::pair(tup1, tup2))] or pair_map[(std::pair(tup2, tup1))])){
                    condition_two = true;
                }

                

                if(!edge(tup1.first, tup2.first, g1->g).second && !edge(tup1.second, tup2.second, g2->g).second and not condition_two){
                    condition_two = true;   
                }
              


                if(condition_two){ 
                    auto e = add_edge(node1, node2, new_mod_product.g);
                    auto pair1 = std::pair<Vertex, Vertex>(node1, node2);
                    auto pair2 = std::pair<Vertex, Vertex>(node2, node1);
                    new_mod_product.bonds[pair1] = "NONE";
                    new_mod_product.bonds[pair2] = "NONE";
                    

                }
                else if(edge(tup1.first, tup2.first, g1->g).second && edge(tup1.second, tup2.second, g2->g).second && g1->bonds[std::pair<Vertex, Vertex>(tup1.first, tup2.first)] == g2->bonds[std::pair<Vertex, Vertex>(tup1.second, tup2.second)]){
                    auto e = add_edge(node1, node2, new_mod_product.g);
                    auto pair1 = std::pair<Vertex, Vertex>(node1, node2);
                    auto pair2 = std::pair<Vertex, Vertex>(node2, node1);
                    new_mod_product.bonds[pair1] = g1->bonds[std::pair<Vertex, Vertex>(tup1.first, tup2.first)];
                    new_mod_product.bonds[pair2] = g1->bonds[std::pair<Vertex, Vertex>(tup1.first, tup2.first)];
                }




            }

        }


        results.push_back(new_mod_product);



        
    }
    
    std::sort(results.begin(), results.end(),
    [](const ModularProductGraph& a, const ModularProductGraph& b) {
        return num_vertices(a.g) > num_vertices(b.g); 
        }
    );


    return results;


} 


//Finds the type-A connected components of a modular product graph, WITH the removal of certain type-0 edges
std::vector<struct ModularProductGraph> alternate_find_type1_connected_components(struct ModularProductGraph& graph){


    
    
    std::vector<struct ModularProductGraph> results;


    std::vector<int> component (boost::num_vertices (graph.g));
    int num_components = boost::connected_components (graph.g, &component[0]);

    

    std::map<int, std::vector<Vertex>> component_map = std::map<int, std::vector<Vertex>>();


    for(int i=0; i<num_components; i++){
        component_map[i] = std::vector<Vertex>();
    }

    for (const auto& vi : boost::make_iterator_range(boost::vertices(graph.g))){
        component_map[component[vi]].push_back(vi);

    }

    auto g1 = graph.g1;
    auto g2 = graph.g2;
    

    std::map<Vertex, Vertex> inverse_clique_map;

    for (const auto& vi : boost::make_iterator_range(boost::vertices(g1->g))){
        inverse_clique_map[g1->clique_map[vi]] = vi;
    }

    int max_comp = 0;

    for(int i=0; i<num_components; i++){
        if(component_map[i].size()>max_comp){
            max_comp = component_map[i].size();
        }
    }

   

   
    for(int i=0; i<num_components; i++){


        auto curr_vertices = component_map[i];

        if(curr_vertices.size()<maximum){
            continue;
        }

        struct ModularProductGraph new_mod_product;

        new_mod_product.g1 = graph.g1;
        new_mod_product.g2 = graph.g2;

        std::vector<Vertex> good_vertices1;
        std::vector<Vertex> good_vertices2;

        std::map<std::pair<Vertex,Vertex>, Vertex> node_map_to_mod_node;

        std::map<std::pair<Vertex, Vertex>, Vertex> pair_vertex_map;

        for(auto vertex: curr_vertices){
            
            auto v = add_vertex(new_mod_product.g);
            new_mod_product.node_map[v] = graph.node_map[vertex];
            node_map_to_mod_node[graph.node_map[vertex]] = v;
            new_mod_product.atoms[v] = graph.atoms[vertex];
            new_mod_product.clique_map[v] = graph.clique_map[vertex];
            pair_vertex_map[graph.node_map[vertex]] = v;
            if(!is_in(good_vertices1, new_mod_product.node_map[v].first)){
                good_vertices1.push_back(new_mod_product.node_map[v].first);
            }
            if(!is_in(good_vertices2, new_mod_product.node_map[v].second)){
                good_vertices2.push_back(new_mod_product.node_map[v].second);
            }
            
            
        }

        

        std::sort(good_vertices1.begin(), good_vertices1.end());
        std::sort(good_vertices2.begin(), good_vertices2.end());

        

        double num_iter = (num_vertices(new_mod_product.g) * num_vertices(new_mod_product.g))/2;

        int curr_num = 0;
        

        std::set<std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>> pairs;

        std::map<std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>, bool> pair_map;

        int iter = 0;


        
        int it = 0;
        for (const auto& node1 : boost::make_iterator_range(boost::vertices(new_mod_product.g))){
           
            auto tup1 = new_mod_product.node_map[node1];

            int max_num_edges = 0;

            for (const auto& v : boost::make_iterator_range(boost::vertices(new_mod_product.g))){
                if(v<=node1){
                    continue;
                }
                auto tup2 = new_mod_product.node_map[v];

                if(!edge(tup1.first, tup2.first, graph.g1->g).second && !edge(tup1.second, tup2.second, graph.g2->g).second){
                    max_num_edges++;   
                }

            }

            auto feasible_pairs = reachable_pairs(tup1.first, *g1, tup1.second, *g2, good_vertices1, good_vertices2, max_num_edges);

            std::sort(feasible_pairs.begin(), feasible_pairs.end());
             auto new_end = std::unique(feasible_pairs.begin(), feasible_pairs.end());

   
             feasible_pairs.erase(new_end, feasible_pairs.end());
             it++;

            for (auto pair : feasible_pairs){

                auto node2 = node_map_to_mod_node[pair];



                if(node1<=node2){
                    continue;
                }

                iter++;
               
                
                if(edge(node1, node2, new_mod_product.g).second){
                    continue;
                }
                curr_num ++;
                
                
                auto tup2 = new_mod_product.node_map[node2];

                auto edge1 = g1->edge_map[tup1.first];
                auto edge2 = g1->edge_map[tup2.first];
                
                
                if(tup1.first==tup2.first || tup1.second==tup2.second){
                    continue;
                }

                

                auto pair_pair1 = std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>(std::pair<Vertex, Vertex>(graph.g1->clique_map[tup1.first], tup1.second), std::pair<Vertex, Vertex>(graph.g1->clique_map[tup2.first], tup2.second));
                auto pair_pair2 = std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex, Vertex>>(std::pair<Vertex, Vertex>(graph.g1->clique_map[tup2.first], tup2.second), std::pair<Vertex, Vertex>(graph.g1->clique_map[tup1.first], tup1.second));

                
                std::vector<std::pair<std::pair<Vertex, Vertex>, std::pair<Vertex,Vertex>>> new_edges;



                if(!edge(tup1.first, tup2.first, g1->g).second && !edge(tup1.second, tup2.second, g2->g).second){ //!edge(tup1.first, tup2.first, g1->g).second && !edge(tup1.second, tup2.second, g2->g).second and equivalent_paths(tup1.first, tup2.first, *g1, tup1.second, tup2.second, *g2, vertex_pairs)){//(map[pair_pair1])){//(!edge(tup1.first, tup2.first, g1->g).second && !edge(tup1.second, tup2.second, g2->g).second && equivalent_paths(tup1.first, tup2.first, *g1, tup1.second, tup2.second, *g2)){
                    auto e = add_edge(node1, node2, new_mod_product.g);
                    auto pair1 = std::pair<Vertex, Vertex>(node1, node2);
                    auto pair2 = std::pair<Vertex, Vertex>(node2, node1);
                    new_mod_product.bonds[pair1] = "NONE";
                    new_mod_product.bonds[pair2] = "NONE";
                    

                }
                else if(edge(tup1.first, tup2.first, g1->g).second && edge(tup1.second, tup2.second, g2->g).second && g1->bonds[std::pair<Vertex, Vertex>(tup1.first, tup2.first)] == g2->bonds[std::pair<Vertex, Vertex>(tup1.second, tup2.second)]){
                    auto e = add_edge(node1, node2, new_mod_product.g);
                    auto pair1 = std::pair<Vertex, Vertex>(node1, node2);
                    auto pair2 = std::pair<Vertex, Vertex>(node2, node1);
                    new_mod_product.bonds[pair1] = g1->bonds[std::pair<Vertex, Vertex>(tup1.first, tup2.first)];
                    new_mod_product.bonds[pair2] = g1->bonds[std::pair<Vertex, Vertex>(tup1.first, tup2.first)];
                }

            }

        }


        results.push_back(new_mod_product);



        
    }
    
    std::sort(results.begin(), results.end(),
    [](const ModularProductGraph& a, const ModularProductGraph& b) {
        return num_vertices(a.g) > num_vertices(b.g); 
        }
    );


    return results;


} 















template <typename T1, typename T2, typename T3>
struct Triple {
    T1 first;
    T2 second;
    T3 third;

    Triple(T1 f, T2 s, T3 t) : first(f), second(s), third(t) {}
};



//Used for always recursing on largest current clique
bool compareBySetSize(const Triple<std::set<Vertex>, std::set<Vertex>, std::set<Vertex>>& a, const Triple<std::set<Vertex>, std::set<Vertex>, std::set<Vertex>>& b) {
    // Compare based on the size of the first set in each vector
    //return a.third.size() >b.third.size();
    //return a.third.size() + a.second.size() > b.third.size() + b.second.size();
    return a.first.size() < b.first.size();
    
    
}


//Used for always yielding current largest clique
bool compareBySetSize2(const std::set<Vertex>& a, const std::set<Vertex>& b) {
    // Compare based on the size of the first set in each vector
    
    return a.size() < b.size() ;
    
    
}










//Finds all type-A connected cliques in a modular product graph
CliqueGenerator bron_kerbosch_non_recursive(std::set<Vertex> init_R, std::set<Vertex> init_P, std::set<Vertex> init_X, ModularProductGraph &graph, std::vector<std::set<Vertex>> &cliques){

   

    std::list<Triple<std::set<Vertex>, std::set<Vertex>, std::set<Vertex>>> stack;
    

    auto insert_sorted = [](std::list<Triple<std::set<Vertex>, std::set<Vertex>, std::set<Vertex>>>& vec, Triple<std::set<Vertex>, std::set<Vertex>, std::set<Vertex>>& value) {
        auto pos = std::lower_bound(vec.begin(), vec.end(), value,  compareBySetSize); // Find the correct position
        vec.insert(pos, value); // Insert at the position
    };



    auto insert_sorted2 = [](std::vector<std::set<Vertex>>& vec, std::set<Vertex> value) {
        auto pos = std::lower_bound(vec.begin(), vec.end(), value,  compareBySetSize2); // Find the correct position
        vec.insert(pos, value); // Insert at the position
    };

    int cliques_found =0;

    Triple<std::set<Vertex>, std::set<Vertex>, std::set<Vertex>> init(init_R, init_P, init_X);
    

    stack.push_back(init);

    

    int it = 0;



    auto curr = init;


    std::vector<std::set<Vertex>> curr_cliques;

    

    while(stack.size()>0){
        it++;
 
        if(curr_cliques.size()>100){
            
            
            int s = curr_cliques.size();
            for(int i=0; i<s; i++){
                auto cl = curr_cliques.back();
                if(cl.size()>=maximum){
                    co_yield cl;
                }
                curr_cliques.pop_back();
            }
        }
        
        
        curr = stack.back();
        
        
    

        



        

        auto R = curr.first;
        auto P = curr.second;
        auto X = curr.third;

        stack.pop_back();
        
        if(P.size() + R.size() < maximum){
            continue;
        }




        if(P.size()==0 && X.size() == 0){
            cliques_found++;
            if(R.size()>=maximum){
                //co_yield R;
                insert_sorted2(curr_cliques, R);
                continue;
            }
            else{
                continue;
            }
            
        }

        
    
        

        auto ref_P = P;
        
        while(P.size()>0){

            Vertex candidate = -1;

            bool found = false;

            for(auto &pv : P){
                for(auto &rv : R){
                  if(edge(pv, rv, graph.g).second){ 
                    if(graph.bonds[pair(rv,pv)]!="NONE"){
                        candidate = pv;
                        found = true;
                        break;
                    }
                  } 
                }
                if(found){
                    break;
                }
            }

            if(!found ){
                if(R.size()!=0){
                    insert_sorted2(curr_cliques, R);
                    break;
                    
                }
                else{
                    candidate = *P.begin();

                }
                
            }

            auto v = candidate;
           
            std::set<Vertex> new_R = R;
            new_R.insert(v);
            std::set<Vertex> new_P;
            std::set<Vertex> new_X;
            if(R.size()!=0){


                for(auto &vd : X){
                    if(edge(vd, v, graph.g).second){
                        new_X.insert(vd);
                    }
                }
                for(auto &vd : P){
                    if(edge(vd, v, graph.g).second){
                        new_P.insert(vd);
                    }
                }
            }
            else{
                for(auto &vd : X){
                    if(edge(vd, v, graph.g).second){
                        new_X.insert(vd);
                    }
                }
                for(auto &vd : ref_P){
                    if(edge(vd, v, graph.g).second){
                        new_P.insert(vd);
                    }
                }
            }
            struct Triple<std::set<Vertex>, std::set<Vertex>, std::set<Vertex>> new_instance(new_R, new_P, new_X);
            
            
            insert_sorted(stack, new_instance);
            
            P.erase(v);
            X.insert(v);


        }   
           



    }
    int s = curr_cliques.size();
    for(int i=0; i<s; i++){
        auto cl = curr_cliques.back();
        if(cl.size()>=maximum){
            co_yield cl;
        }
        curr_cliques.pop_back();
    }

}




//Driver function that returns a generator that yields all maximal typer-A connected cliques in a modular product graph
CliqueGenerator bron_kerbosch_driver(ModularProductGraph &graph){



    
    std::vector<std::set<Vertex>> cliques;
    std::set<Vertex> P;
    std::set<Vertex> X;
    std::set<Vertex> R;

    for (const auto& node : boost::make_iterator_range(boost::vertices(graph.g))){
        P.insert(node);
    }

    auto clique_generator = bron_kerbosch_non_recursive(R,P,X, graph, cliques);

    return clique_generator;



}




//Parses a smiles string and constructs a simple, labeled graph
MolecularGraph smiles_to_graph(std::string smiles){



    struct MolecularGraph g;

    OpenBabel::OBMol mol;
         

    OpenBabel::OBConversion conv;
    

    conv.SetInFormat("smi");

   
    
    conv.ReadString(&mol, smiles);    

    //mol.AddHydrogens();

    auto vertex_map = std::map<int,Vertex>();


    int map_num = 0;

    

     for (OpenBabel::OBAtomIterator atom = mol.BeginAtoms(); atom != mol.EndAtoms(); ++atom){

        auto a = *atom;

        auto id =  a->GetIdx();
        
        Vertex v = add_vertex(g.g);
        vertex_map[id] = v;
        g.atoms[v] = std::to_string(a->GetAtomicNum()) + std::to_string(a->GetFormalCharge());
        

    }


    for (OpenBabel::OBBondIterator bond = mol.BeginBonds(); bond != mol.EndBonds(); ++bond){

        auto b = *bond;

        auto src =  vertex_map[b->GetBeginAtomIdx()];
        auto tar = vertex_map[b->GetEndAtomIdx()];
       auto e = add_edge(src, tar, g.g).first;

      g.bonds[std::pair(src, tar)] = std::to_string(b->GetBondOrder());
      g.bonds[std::pair(tar, src)] = std::to_string(b->GetBondOrder());
      
     


      g.clique_map[std::pair(tar, src)] = map_num;
      g.clique_map[std::pair(src, tar)] = map_num;

      map_num++;
      
      
    }

     



    return g;


}

struct GMLRule{
    // The left and right graphs include the rule context.
    std::string rule_id;
    MolecularGraph left;
    MolecularGraph context;
    MolecularGraph right;
};

struct GMLToken{
    std::string text;
    int line;
    int column;
};

struct GMLValue{
    bool is_list = false;
    std::string scalar;
    std::vector<std::pair<std::string, GMLValue>> entries;
};

struct GMLEdgeRecord{
    int source;
    int target;
    std::string label;
};

struct GMLGraphSpec{
    std::map<int, std::string> nodes;
    std::vector<GMLEdgeRecord> edges;
};

static std::string gml_trim(std::string value){
    auto first = std::find_if(value.begin(), value.end(), [](unsigned char c){
        return !std::isspace(c);
    });
    auto last = std::find_if(value.rbegin(), value.rend(), [](unsigned char c){
        return !std::isspace(c);
    }).base();

    if(first >= last){
        return "";
    }
    return std::string(first, last);
}

static std::string gml_lower(std::string value){
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c){
        return std::tolower(c);
    });
    return value;
}

static bool gml_key_equals(const std::string& a, const std::string& b){
    return gml_lower(a) == gml_lower(b);
}

static std::vector<GMLToken> tokenize_gml(const std::string& contents){
    std::vector<GMLToken> tokens;
    int line = 1;
    int column = 1;
    std::size_t i = 0;

    auto advance = [&](){
        if(contents[i] == '\n'){
            line++;
            column = 1;
        }
        else{
            column++;
        }
        i++;
    };

    while(i < contents.size()){
        char current = contents[i];

        if(std::isspace(static_cast<unsigned char>(current))){
            advance();
            continue;
        }

        if(current == '#'){
            while(i < contents.size() && contents[i] != '\n'){
                advance();
            }
            continue;
        }

        if(current == '[' || current == ']'){
            tokens.push_back({std::string(1, current), line, column});
            advance();
            continue;
        }

        if(current == '"'){
            int start_line = line;
            int start_column = column;
            std::string value;
            advance();

            bool closed = false;
            while(i < contents.size()){
                current = contents[i];
                if(current == '"'){
                    advance();
                    closed = true;
                    break;
                }

                if(current == '\\'){
                    advance();
                    if(i >= contents.size()){
                        break;
                    }

                    char escaped = contents[i];
                    if(escaped == 'n'){
                        value.push_back('\n');
                    }
                    else if(escaped == 't'){
                        value.push_back('\t');
                    }
                    else{
                        value.push_back(escaped);
                    }
                    advance();
                    continue;
                }

                value.push_back(current);
                advance();
            }

            if(!closed){
                throw std::runtime_error("Unterminated quoted string in GML file at line " + std::to_string(start_line) + ", column " + std::to_string(start_column));
            }

            tokens.push_back({value, start_line, start_column});
            continue;
        }

        int start_line = line;
        int start_column = column;
        std::string value;
        while(i < contents.size()){
            current = contents[i];
            if(std::isspace(static_cast<unsigned char>(current)) || current == '[' || current == ']' || current == '"' || current == '#'){
                break;
            }
            value.push_back(current);
            advance();
        }

        if(!value.empty()){
            tokens.push_back({value, start_line, start_column});
        }
    }

    return tokens;
}

class GMLParser{
public:
    GMLParser(std::vector<GMLToken> tokens) : tokens(std::move(tokens)) {}

    GMLValue parse_document(){
        GMLValue document;
        document.is_list = true;

        while(!eof()){
            document.entries.push_back(parse_entry());
        }

        return document;
    }

private:
    std::vector<GMLToken> tokens;
    std::size_t position = 0;

    bool eof() const{
        return position >= tokens.size();
    }

    const GMLToken& peek() const{
        return tokens[position];
    }

    const GMLToken& consume(){
        return tokens[position++];
    }

    bool accept(const std::string& text){
        if(!eof() && peek().text == text){
            position++;
            return true;
        }
        return false;
    }

    std::pair<std::string, GMLValue> parse_entry(){
        if(eof()){
            throw std::runtime_error("Unexpected end of GML file while reading a key");
        }

        const auto key = consume();
        if(key.text == "[" || key.text == "]"){
            throw std::runtime_error("Expected a GML key at line " + std::to_string(key.line) + ", column " + std::to_string(key.column));
        }

        return {key.text, parse_value()};
    }

    GMLValue parse_value(){
        if(eof()){
            throw std::runtime_error("Unexpected end of GML file while reading a value");
        }

        if(accept("[")){
            GMLValue value;
            value.is_list = true;

            while(!accept("]")){
                if(eof()){
                    throw std::runtime_error("Unterminated GML list");
                }
                value.entries.push_back(parse_entry());
            }

            return value;
        }

        const auto value_token = consume();
        if(value_token.text == "]"){
            throw std::runtime_error("Unexpected ']' at line " + std::to_string(value_token.line) + ", column " + std::to_string(value_token.column));
        }

        GMLValue value;
        value.scalar = value_token.text;
        return value;
    }
};

static GMLValue parse_gml_document(const std::string& filename){
    std::ifstream input(filename);
    if(!input){
        throw std::runtime_error("Could not open GML file: " + filename);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();

    GMLParser parser(tokenize_gml(buffer.str()));
    return parser.parse_document();
}

static const GMLValue* gml_find_first(const GMLValue& value, const std::string& key){
    if(!value.is_list){
        return nullptr;
    }

    for(const auto& entry : value.entries){
        if(gml_key_equals(entry.first, key)){
            return &entry.second;
        }
    }

    return nullptr;
}

static std::vector<const GMLValue*> gml_find_all(const GMLValue& value, const std::string& key){
    std::vector<const GMLValue*> matches;
    if(!value.is_list){
        return matches;
    }

    for(const auto& entry : value.entries){
        if(gml_key_equals(entry.first, key)){
            matches.push_back(&entry.second);
        }
    }

    return matches;
}

static bool gml_get_scalar(const GMLValue& value, std::initializer_list<std::string> keys, std::string& result){
    for(const auto& key : keys){
        auto child = gml_find_first(value, key);
        if(child == nullptr){
            continue;
        }
        if(child->is_list){
            throw std::runtime_error("Expected scalar GML value for key '" + key + "'");
        }
        result = child->scalar;
        return true;
    }

    return false;
}

static int gml_parse_int(const std::string& value, const std::string& context){
    std::string trimmed = gml_trim(value);
    try{
        std::size_t parsed = 0;
        int result = std::stoi(trimmed, &parsed);
        if(parsed != trimmed.size()){
            throw std::runtime_error("");
        }
        return result;
    }
    catch(...){
        throw std::runtime_error("Expected integer for " + context + ", got '" + value + "'");
    }
}

static int gml_required_int(const GMLValue& value, std::initializer_list<std::string> keys, const std::string& context){
    std::string result;
    if(!gml_get_scalar(value, keys, result)){
        throw std::runtime_error("Missing integer field in " + context);
    }
    return gml_parse_int(result, context);
}

static int gml_atomic_number(const std::string& symbol){
    static const std::map<std::string, int> atomic_numbers = {
        {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8},
        {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16},
        {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24},
        {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32},
        {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40},
        {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48},
        {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55}, {"Ba", 56},
        {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60}, {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64},
        {"Tb", 65}, {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, {"Lu", 71}, {"Hf", 72},
        {"Ta", 73}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80},
        {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88},
        {"Ac", 89}, {"Th", 90}, {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95}, {"Cm", 96},
        {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100}, {"Md", 101}, {"No", 102}, {"Lr", 103},
        {"Rf", 104}, {"Db", 105}, {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110},
        {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114}, {"Mc", 115}, {"Lv", 116}, {"Ts", 117}, {"Og", 118}
    };

    auto found = atomic_numbers.find(symbol);
    if(found == atomic_numbers.end()){
        return 0;
    }
    return found->second;
}

static std::string gml_normalize_symbol(std::string symbol){
    if(symbol.empty()){
        return symbol;
    }

    symbol[0] = std::toupper(static_cast<unsigned char>(symbol[0]));
    for(std::size_t i = 1; i < symbol.size(); i++){
        symbol[i] = std::tolower(static_cast<unsigned char>(symbol[i]));
    }
    return symbol;
}

static bool gml_parse_charge_suffix(const std::string& suffix, int& charge){
    if(suffix.empty()){
        charge = 0;
        return true;
    }

    if(suffix[0] != '+' && suffix[0] != '-'){
        return false;
    }

    int sign = suffix[0] == '+' ? 1 : -1;
    bool repeated_signs = true;
    for(char c : suffix){
        if(c != suffix[0]){
            repeated_signs = false;
            break;
        }
    }

    if(repeated_signs){
        charge = sign * static_cast<int>(suffix.size());
        return true;
    }

    std::string magnitude = suffix.substr(1);
    if(magnitude.empty()){
        charge = sign;
        return true;
    }

    for(char c : magnitude){
        if(!std::isdigit(static_cast<unsigned char>(c))){
            return false;
        }
    }

    charge = sign * std::stoi(magnitude);
    return true;
}

static std::string normalize_gml_atom_label(std::string label){
    label = gml_trim(label);
    if(label.size() >= 2 && label.front() == '[' && label.back() == ']'){
        label = gml_trim(label.substr(1, label.size() - 2));
    }

    if(label.empty()){
        return label;
    }

    for(int symbol_length = std::min<int>(2, label.size()); symbol_length >= 1; symbol_length--){
        std::string symbol = gml_normalize_symbol(label.substr(0, symbol_length));
        int atomic_number = gml_atomic_number(symbol);
        if(atomic_number == 0){
            continue;
        }

        int charge = 0;
        if(gml_parse_charge_suffix(label.substr(symbol_length), charge)){
            return std::to_string(atomic_number) + std::to_string(charge);
        }
    }

    return label;
}

static std::string normalize_gml_bond_label(std::string label){
    label = gml_trim(label);
    std::string lower = gml_lower(label);

    if(lower == "-" || lower == "single"){
        return "1";
    }
    if(lower == "=" || lower == "double"){
        return "2";
    }
    if(lower == "#" || lower == "triple"){
        return "3";
    }
    if(lower == ":" || lower == "aromatic"){
        return "5";
    }

    return label;
}

static std::string gml_node_label(const GMLValue& node){
    std::string atomic_number;
    if(gml_get_scalar(node, {"atomicNumber", "atomic_number", "atomicnum"}, atomic_number)){
        std::string charge_string;
        int charge = 0;
        if(gml_get_scalar(node, {"charge", "formalCharge", "formal_charge"}, charge_string)){
            charge = gml_parse_int(charge_string, "node charge");
        }
        return std::to_string(gml_parse_int(atomic_number, "node atomic number")) + std::to_string(charge);
    }

    std::string label;
    if(gml_get_scalar(node, {"label", "atom", "element", "type"}, label)){
        return normalize_gml_atom_label(label);
    }

    return "";
}

static std::string gml_edge_label(const GMLValue& edge){
    std::string label;
    if(gml_get_scalar(edge, {"label", "bond", "order", "type"}, label)){
        return normalize_gml_bond_label(label);
    }
    return "";
}

static const GMLValue& gml_graph_body(const GMLValue& value){
    if(gml_find_first(value, "node") != nullptr || gml_find_first(value, "edge") != nullptr){
        return value;
    }

    auto graph = gml_find_first(value, "graph");
    if(graph != nullptr && graph->is_list){
        return *graph;
    }

    return value;
}

static GMLGraphSpec parse_gml_graph_spec(const GMLValue& value, const std::string& context){
    const GMLValue& body = gml_graph_body(value);
    if(!body.is_list){
        throw std::runtime_error("Expected GML list for " + context);
    }

    GMLGraphSpec spec;
    for(const auto node : gml_find_all(body, "node")){
        if(!node->is_list){
            throw std::runtime_error("Expected GML node list in " + context);
        }

        int id = gml_required_int(*node, {"id"}, context + " node id");
        std::string label = gml_node_label(*node);
        auto existing = spec.nodes.find(id);
        if(existing != spec.nodes.end() && existing->second != label){
            throw std::runtime_error("Conflicting labels for node " + std::to_string(id) + " in GML " + context);
        }
        spec.nodes[id] = label;
    }

    for(const auto edge : gml_find_all(body, "edge")){
        if(!edge->is_list){
            throw std::runtime_error("Expected GML edge list in " + context);
        }

        GMLEdgeRecord record;
        record.source = gml_required_int(*edge, {"source", "src", "from"}, context + " edge source");
        record.target = gml_required_int(*edge, {"target", "tgt", "to"}, context + " edge target");
        record.label = gml_edge_label(*edge);
        spec.edges.push_back(record);
    }

    return spec;
}

static GMLGraphSpec merge_gml_graph_specs(const GMLGraphSpec& context, const GMLGraphSpec& delta, const std::string& section_name){
    GMLGraphSpec merged = context;

    for(const auto& node : delta.nodes){
        auto existing = merged.nodes.find(node.first);
        if(existing != merged.nodes.end() && !existing->second.empty() && !node.second.empty() && existing->second != node.second){
            throw std::runtime_error("Conflicting labels for node " + std::to_string(node.first) + " while building GML " + section_name + " graph");
        }

        if(existing == merged.nodes.end() || existing->second.empty()){
            merged.nodes[node.first] = node.second;
        }
    }

    for(const auto& edge : delta.edges){
        merged.edges.push_back(edge);
    }

    return merged;
}

static Vertex ensure_gml_vertex(MolecularGraph& graph, std::map<int, Vertex>& vertex_map, const GMLGraphSpec& spec, int id){
    auto existing = vertex_map.find(id);
    if(existing != vertex_map.end()){
        return existing->second;
    }

    Vertex vertex = add_vertex(graph.g);
    vertex_map[id] = vertex;

    auto label = spec.nodes.find(id);
    graph.atoms[vertex] = label == spec.nodes.end() ? "" : label->second;
    return vertex;
}

static MolecularGraph gml_graph_spec_to_molecular_graph(const GMLGraphSpec& spec, const std::string& context){
    MolecularGraph graph;
    std::map<int, Vertex> vertex_map;

    for(const auto& node : spec.nodes){
        ensure_gml_vertex(graph, vertex_map, spec, node.first);
    }

    std::map<std::pair<int, int>, std::string> seen_edges;
    int edge_index = 0;

    for(const auto& edge_record : spec.edges){
        int first_id = std::min(edge_record.source, edge_record.target);
        int second_id = std::max(edge_record.source, edge_record.target);
        std::pair<int, int> id_pair(first_id, second_id);

        auto existing = seen_edges.find(id_pair);
        if(existing != seen_edges.end()){
            if(existing->second != edge_record.label){
                throw std::runtime_error("Conflicting labels for edge " + std::to_string(edge_record.source) + "-" + std::to_string(edge_record.target) + " in GML " + context);
            }
            continue;
        }
        seen_edges[id_pair] = edge_record.label;

        if(spec.nodes.find(edge_record.source) == spec.nodes.end() || spec.nodes.find(edge_record.target) == spec.nodes.end()){
            throw std::runtime_error("GML " + context + " edge " + std::to_string(edge_record.source) + "-" + std::to_string(edge_record.target) + " references an undeclared node");
        }

        Vertex source_vertex = ensure_gml_vertex(graph, vertex_map, spec, edge_record.source);
        Vertex target_vertex = ensure_gml_vertex(graph, vertex_map, spec, edge_record.target);
        add_edge(source_vertex, target_vertex, graph.g);

        graph.bonds[std::pair(source_vertex, target_vertex)] = edge_record.label;
        graph.bonds[std::pair(target_vertex, source_vertex)] = edge_record.label;
        graph.clique_map[std::pair(source_vertex, target_vertex)] = edge_index;
        graph.clique_map[std::pair(target_vertex, source_vertex)] = edge_index;
        edge_index++;
    }

    return graph;
}

MolecularGraph parse_gml_graph_file(std::string gml_file){
    GMLValue document = parse_gml_document(gml_file);
    const GMLValue* graph = gml_find_first(document, "graph");
    const GMLValue& graph_value = graph == nullptr ? document : *graph;
    return gml_graph_spec_to_molecular_graph(parse_gml_graph_spec(graph_value, "graph"), "graph");
}

GMLRule parse_gml_rule_file(std::string gml_file){
    GMLValue document = parse_gml_document(gml_file);
    const GMLValue* rule_value = gml_find_first(document, "rule");
    const GMLValue& rule_body = rule_value == nullptr ? document : *rule_value;

    const GMLValue* left_value = gml_find_first(rule_body, "left");
    const GMLValue* context_value = gml_find_first(rule_body, "context");
    const GMLValue* right_value = gml_find_first(rule_body, "right");
    std::string rule_id;
    gml_get_scalar(rule_body, {"ruleID", "ruleId", "rule_id", "name"}, rule_id);

    if(left_value == nullptr && context_value == nullptr && right_value == nullptr){
        const GMLValue* graph_value = gml_find_first(rule_body, "graph");
        MolecularGraph graph = graph_value == nullptr
            ? parse_gml_graph_file(gml_file)
            : gml_graph_spec_to_molecular_graph(parse_gml_graph_spec(*graph_value, "graph"), "graph");
        return {rule_id, graph, graph, graph};
    }

    GMLGraphSpec context_spec = context_value == nullptr ? GMLGraphSpec() : parse_gml_graph_spec(*context_value, "context");
    GMLGraphSpec left_delta = left_value == nullptr ? GMLGraphSpec() : parse_gml_graph_spec(*left_value, "left");
    GMLGraphSpec right_delta = right_value == nullptr ? GMLGraphSpec() : parse_gml_graph_spec(*right_value, "right");

    GMLRule rule;
    rule.rule_id = rule_id;
    rule.context = gml_graph_spec_to_molecular_graph(context_spec, "context");
    rule.left = gml_graph_spec_to_molecular_graph(merge_gml_graph_specs(context_spec, left_delta, "left"), "left");
    rule.right = gml_graph_spec_to_molecular_graph(merge_gml_graph_specs(context_spec, right_delta, "right"), "right");
    return rule;
}


//Converts a modular product graph of line graphs to its original simple, labeled graph
MolecularGraph get_original_graph_from_modular_linegraph(struct ModularProductGraph& graph){


    struct MolecularGraph orig_graph;


    std::set<Vertex> vertices;

    std::set<std::pair<Vertex, Vertex>> edges;


    for (auto node : boost::make_iterator_range(boost::vertices(graph.g))){
        auto v = graph.node_map[node].first;

        if(vertices.find(graph.g1->edge_map[v].first)==vertices.end()){

            vertices.insert(graph.g1->edge_map[v].first);
        }

        if(vertices.find(graph.g1->edge_map[v].second)==vertices.end()){
            vertices.insert(graph.g1->edge_map[v].second);
        }
        

        edges.insert(graph.g1->edge_map[v]);
    }



    std::map<Vertex, Vertex> vertex_map;

    for(auto vertex : vertices){
        auto v = add_vertex(orig_graph.g);
        vertex_map[vertex] = v;
        orig_graph.atoms[v] = graph.g1->org_graph->atoms[vertex];


    }
    for (auto edge : boost::make_iterator_range(boost::edges(graph.g1->org_graph->g))){
        
    
        if(edges.find(std::pair<Vertex, Vertex>(edge.m_source, edge.m_target))!=edges.end()  or edges.find(std::pair<Vertex, Vertex>(edge.m_target, edge.m_source))!=edges.end()){
           
            auto e = add_edge(vertex_map[edge.m_source], vertex_map[edge.m_target], orig_graph.g);
            orig_graph.bonds[std::pair<Vertex, Vertex>(vertex_map[edge.m_source], vertex_map[edge.m_target])] = graph.g1->org_graph->bonds[(std::pair<Vertex, Vertex>(edge.m_source, edge.m_target))];
            orig_graph.bonds[std::pair<Vertex, Vertex>(vertex_map[edge.m_target], vertex_map[edge.m_source])] = graph.g1->org_graph->bonds[(std::pair<Vertex, Vertex>(edge.m_source, edge.m_target))];

            orig_graph.clique_map[std::pair<Vertex, Vertex>(vertex_map[edge.m_target], vertex_map[edge.m_source])] = graph.g1->org_graph->clique_map[(std::pair<Vertex, Vertex>(edge.m_source, edge.m_target))];
            orig_graph.clique_map[std::pair<Vertex, Vertex>(vertex_map[edge.m_source], vertex_map[edge.m_target])] = graph.g1->org_graph->clique_map[(std::pair<Vertex, Vertex>(edge.m_source, edge.m_target))];
           
        }



    }


    return orig_graph;



}


//Finds the maximum common induced subgraph of a list of line graphs
//Stores results in all_maximum_graphs
void maximum_common_subgraph(std::vector<struct LineGraph> &graphs, bool first){
    
    

    if(first){
        max_size_type1 = 0;
        path_finding_time = 0;
        num_cliques_found = 0;
        maximum = 0;
        maximum_graph = MolecularGraph();
        all_cliques = std::vector<std::vector<Vertex>>();
        num_graphs = graphs.size();

        std::vector<std::vector<std::vector<Vertex>>> outer_vector;

        for(int i=0; i<graphs.size(); i++){
            std::vector<std::vector<Vertex>> inner_vector;
            outer_vector.push_back(inner_vector);
            graphs[i].num = i;
        }
        clique_vectors = outer_vector;
        edges_added = 0;
        edges_removed = 0;
        all_maximum_graphs = std::vector<MolecularGraph>();

        clique_counter = std::vector<int>();

        for(int i=1; i<graphs.size(); i++){
            clique_counter.push_back(0);
        }
        
        

        

        



    }
    if(num_vertices(graphs[0].g)<maximum){
        return;
    }


    if(graphs.size()==2){
        auto mcs_gen = maximal_common_subgraphs(graphs[0], graphs[1]);
        while (mcs_gen.move_next()) {
            MolecularGraph graph = mcs_gen.current_value();

            

           
            if(num_edges(graph.g)==maximum){
                all_maximum_graphs.push_back(graph);
            }
            else if(num_edges(graph.g)>maximum){
                maximum = num_edges(graph.g);
                all_maximum_graphs = std::vector<MolecularGraph>();
                maximum_graph = graph;
                all_maximum_graphs.push_back(maximum_graph);
            }
        }

        

    }
    else{

        auto mcs_gen = maximal_common_subgraphs(graphs[0], graphs[1]);
        while (mcs_gen.move_next()) {

            MolecularGraph subgraph = mcs_gen.current_value();


            if(num_edges(subgraph.g)<maximum){
                continue;
            }

            
            auto sg = molecular_graph_to_line_graph(subgraph);

            std::vector<LineGraph> new_graphs;
            new_graphs.push_back(sg);
            for(int i=2; i<graphs.size(); i++){
                new_graphs.push_back(graphs[i]);
            }
            maximum_common_subgraph(new_graphs, false);
        }
        
       
    }


    

}



//Computes the minmax similarities between a list of line graphs.
//Name of function doesn't make sense
std::map<std::pair<int,int>, int> strong_product_kernel(std::vector<LineGraph> &graphs){

    std::map<std::pair<int,int>, int> kernel;

    double min = 9999999;

    for(int i=0; i<graphs.size(); i++){
        for(int j=0; j<graphs.size(); j++){
            if(i<=j){
                continue;
            }
            auto strong_product = modular_graph_product(graphs[i], graphs[j]);



            std::vector<int> component(boost::num_vertices(strong_product.g));
    
            // Compute connected components
            int num_components = boost::connected_components(strong_product.g, &component[0]);

            // Vector to store the size of each component
            std::vector<std::size_t> component_sizes(num_components, 0);



            // Count vertices in each component
            for(int k = 0; k<component.size(); k++){
                ++component_sizes[component[k]];
               // component_sizes[component[k]] += degree(k, strong_product.g);
            }

            /*
            for (const auto& comp : component) {
                ++component_sizes[comp];

            }*/



            // Find and return the size of the largest connected component
            double max_size = 0; //*std::max_element(component_sizes.begin(), component_sizes.end());

            for(auto s : component_sizes){
                if(s>max_size){
                    max_size = s;
                }
            }
            
            if(max_size<min){
                min = max_size;
            }

            kernel[std::pair(i,j)] = max_size;
            kernel[std::pair(j,i)] = max_size;

        }
    }

    
    return kernel;


}



//Takes a list of line graphs and orders them greedily according to minmax similarity
Triple<std::vector<LineGraph>, std::vector<int>, int> find_ordering(std::vector<LineGraph> &graphs){

    std::vector<LineGraph> ordering;

    std::vector<int> index_ordering;

    std::pair<int,int> min_coords;

    

    auto kernel = strong_product_kernel(graphs);
    
    int min = 9999999;
    for(int i=0; i<graphs.size(); i++){
        for(int j=0; j<graphs.size(); j++){
            if(i<=j){
                continue;
            }

            if(kernel[std::pair(i,j)]<min){
                min = kernel[std::pair(i,j)];
                min_coords.first = i;
                min_coords.second = j;
            }

        }
    }

    int final_min = min;
   
    index_ordering.push_back(min_coords.first);
    index_ordering.push_back(min_coords.second);

    while(index_ordering.size()<graphs.size()){
        min = 9999999;
        int min_coord = 0;
        for(int i=0; i<graphs.size(); i++){
            if(is_in(index_ordering, i)){
                continue;
            }
            int max = 0;
            for(int j : index_ordering){
                if(kernel[std::pair(i,j)]>max){
                    max = kernel[std::pair(i,j)];
                }
            }
            if(max<min){
                min = max;
                min_coord = i;
            }
            

        }
        index_ordering.push_back(min_coord);


    }

    for(int i : index_ordering){
        ordering.push_back(graphs[i]);
    }







    return Triple(ordering, index_ordering, final_min);




}


// Splits a string according to delimiter
std::vector<std::string> split(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}





//Takes a molecular graph and returns a corresponding smiles string
std::string to_mol(MolecularGraph &g){

    OpenBabel::OBMol mol;

    std::map<Vertex, int> atom_map;

    int index = 1;

    for (auto vertex : boost::make_iterator_range(boost::vertices(g.g))){
        std::string type = g.atoms[vertex];
        int charge = type.back() - '0';
        type.pop_back();
        OpenBabel::OBAtom* c1 = mol.NewAtom();
        c1->SetAtomicNum(std::stoi(type)); 
        if(charge!=0){
            c1->SetFormalCharge(charge);
        }
        
        atom_map[vertex] = index;
        index ++;

    }
    for (auto edge : boost::make_iterator_range(boost::edges(g.g))){

        std::string bond_type = g.bonds[(std::pair(source(edge, g.g), target(edge, g.g)))];



        int int_bond = std::stoi(bond_type);

        mol.AddBond(atom_map[source(edge, g.g)], atom_map[target(edge, g.g)], int_bond);




    }

    mol.ConnectTheDots();
    //mol.PerceiveBondOrders();

    // Output as SMILES
    OpenBabel::OBConversion conv;
    conv.SetOutFormat("smi");
    std::string smiles = conv.WriteString(&mol);
    
    smiles.pop_back();
    smiles.pop_back();

    return smiles;




}


//Takes input csv file and solves all instances, outputting results in output csv file
void solve_instances(std::string input_file, std::string output_file, bool prune, bool minmax_order){

    std::vector<std::vector<std::string>> instances;

     

     ofstream results(output_file);

     

    std::ifstream infile(input_file);
    std::string line;
    while (std::getline(infile, line))
    {

        std::istringstream iss(line);

        std::vector<std::string> instance;

        for(auto s : split(line, ",")){
            if(s.size()>1){
                instance.push_back(s);
            }    
        }
        instances.push_back(instance);

         
          
        
    }

    int iter = 1;

    for(auto instance : instances){


        std::vector<LineGraph> graphs;
        std::vector<MolecularGraph> mgraphs;

        
        for(auto smiles : instance){           
           

            auto g = smiles_to_graph(smiles);
            
             
            auto lgg =  molecular_graph_to_line_graph(g);

            graphs.push_back(lgg);
            

        }

        if(minmax_order){
            auto ordering = find_ordering(graphs); 
        
            graphs = ordering.first;
        }

        

        special_product = prune;

        maximum_common_subgraph(graphs, true);

        for(auto mcs_graph : all_maximum_graphs){
            results << to_mol(mcs_graph) << ",";
        }
        results << endl;

        cout << "Instance: " << iter << "/" << instances.size() << " solved." << endl;

        iter++;






    }

}



/**********************functions for experiments*************************************************************************************/




void find_all_cliques_per_file(){

    std::ifstream infile("smiles_list9");
    std::vector<std::string> smiles;
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        smiles.push_back(line);       
        
    } 

    ofstream MyFile("zinc12k_common_subgraph_data_smiles_list9");
    int iter = 0;
    for(auto data_line : smiles){


            auto split_string = split(data_line, ",");
            

            int k = stoi(split_string[0]);
            int j = stoi(split_string[1]);
            auto sm1 = split_string[2];
            auto sm2 = split_string[3];
            if(k<j){
                continue;
            }  
            if(iter<502){
                //iter++;
                //continue;
            } 
            iter++;

            std::map<int, int> clique_count;

            for(int i=1; i<11; i++){
                clique_count[i] = 0;
            }

            auto g1 = smiles_to_graph(sm1);
            auto lg1 = molecular_graph_to_line_graph(g1);
            
            auto g2 = smiles_to_graph(sm2);
            auto lg2 = molecular_graph_to_line_graph(g2);


            auto direct_product = modular_graph_product(lg1, lg2);

            special_product = false;

            auto type1_connected_components = alternate_find_type1_connected_components(direct_product);

            cout << iter << ": DONE!" << endl;

            auto now = chrono::system_clock::now();

            auto duration = now.time_since_epoch();

            auto milliseconds1
                = chrono::duration_cast<chrono::nanoseconds>(
                      duration)
                      .count();

            int cliques_num = 0;

            for(auto subgraph : type1_connected_components){


                

                auto cliques =  bron_kerbosch_driver(subgraph);
                //for(const auto& clique: cliques){
                while (cliques.move_next()){

                    auto clique = cliques.current_value();

                    clique_count[clique.size()]++;
                }

                

            

            }

            

            MyFile << k << "," <<sm1 << endl;
            MyFile << j << "," <<sm2 << endl;

            for(int i=1; i<11; i++){
                MyFile << i << ":" << clique_count[i] << endl;
            }
            MyFile << "#" << endl;

        
    }







}


void find_all_cliques(){

    std::ifstream infile("smiles_mcs.txt");
    std::vector<std::string> smiles;
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        smiles.push_back(line);       
        
    } 

    ofstream MyFile("zinc250k_common_subgraph_data.txt");
    int iter = 0;
    for(int k=0; k<500; k++){
        for(int j=0; j<500; j++){
            auto sm1 = smiles[k];
            auto sm2 = smiles[j];
            if(k<j){
                continue;
            }  
            
            iter++;

            std::map<int, int> clique_count;

            for(int i=1; i<100; i++){
                clique_count[i] = 0;
            }

            auto g1 = smiles_to_graph(sm1);
            auto lg1 = molecular_graph_to_line_graph(g1);
            
            auto g2 = smiles_to_graph(sm2);
            auto lg2 = molecular_graph_to_line_graph(g2);


            auto direct_product = modular_graph_product(lg1, lg2);

            special_product = true;

            auto type1_connected_components = alternate_find_type1_connected_components(direct_product);

            cout << iter << ": DONE!" << endl;

            auto now = chrono::system_clock::now();

            auto duration = now.time_since_epoch();

            auto milliseconds1
                = chrono::duration_cast<chrono::nanoseconds>(
                      duration)
                      .count();

            int cliques_num = 0;

            for(auto subgraph : type1_connected_components){


                

                auto cliques =  bron_kerbosch_driver(subgraph);
                //for(const auto& clique: cliques){
                while (cliques.move_next()){

                    auto clique = cliques.current_value();

                    clique_count[clique.size()]++;
                }

                

            

            }

            

            MyFile << k << "," <<sm1 << endl;
            MyFile << j << "," <<sm2 << endl;

            for(int i=1; i<100; i++){
                MyFile << i << ":" << clique_count[i] << endl;
            }
            MyFile << "#" << endl;

        }
    }






}


void test_cliques_passed(){
    std::ifstream infile("test_smiles3");
    std::vector<std::string> smiles;
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        smiles.push_back(line);       
        
    }
    std::vector<std::string> smiles2;
    std::ifstream infile2("test_smiles27");
    while (std::getline(infile2, line))
    {
        std::istringstream iss(line);
        smiles2.push_back(line);       
        
    }
    std::mt19937_64 random_engine{std::random_device{}()};

    int iter = 0;
    ofstream MyFile("special_product_results.txt");
    while(true){
        std::vector<std::string> instance{};
        std::ranges::sample(smiles, std::back_inserter(instance), 4, random_engine);

        
        std::vector<LineGraph> graphs;
        auto g = smiles_to_graph(instance[0]);
        auto lg = molecular_graph_to_line_graph(g);
        graphs.push_back(lg);

        auto h = smiles_to_graph(instance[1]);
        auto hg = molecular_graph_to_line_graph(h);
        graphs.push_back(hg);
        
        auto q = smiles_to_graph(instance[2]);
        auto  qg = molecular_graph_to_line_graph(q);
        graphs.push_back(qg);
        
        

        std::vector<std::string> instance2{};
        std::ranges::sample(smiles2, std::back_inserter(instance2), 3, random_engine);

        auto w = smiles_to_graph(instance2[0]);
        auto  wg = molecular_graph_to_line_graph(w);
        graphs.push_back(wg);

        auto r = smiles_to_graph(instance2[1]);
        auto  rg = molecular_graph_to_line_graph(r);
        graphs.push_back(rg);

        auto rng = std::default_random_engine {};
        std::shuffle(std::begin(graphs), std::end(graphs), rng);

         auto now = chrono::system_clock::now();

    auto duration = now.time_since_epoch();

    auto milliseconds1
        = chrono::duration_cast<chrono::nanoseconds>(
              duration)
              .count();


        MyFile << "SPECIAL_" << iter << ":" << endl;
        
        special_product = false;

        maximum_common_subgraph(graphs, true);

        MyFile << "    "<<"MCS SIZE: " <<  num_edges(maximum_graph.g) << endl;

        MyFile << "    " <<"CLIQUES PASSED: ";

        for(auto i : clique_counter){
            MyFile << i << ", ";
        }
        MyFile << endl;

        MyFile << "    " <<"CLIQUES FOUND: " << num_cliques_found << endl;


        now = chrono::system_clock::now();
       // Convert the current time to time since epoch
        auto duration2 = now.time_since_epoch();

        auto milliseconds2
            = chrono::duration_cast<chrono::nanoseconds>(
                  duration2)
                  .count();

        auto milliseconds = (double(milliseconds2)- double(milliseconds1))/1000000000;

       MyFile << "    " <<"TIME: " << milliseconds << endl;

       MyFile << "    " << "PATH FINDING TIME: " << path_finding_time << endl;





      now = chrono::system_clock::now();

    duration = now.time_since_epoch();

    milliseconds1
        = chrono::duration_cast<chrono::nanoseconds>(
              duration)
              .count();



        special_product = false;

        MyFile << "NORMAL_" << iter  << ":" << endl;

        maximum_common_subgraph(graphs, true);

        MyFile << "    " <<"MCS SIZE: " << num_edges(maximum_graph.g) << endl;

        MyFile << "    " <<"CLIQUES PASSED: ";

        for(auto i : clique_counter){
            MyFile << i << ", ";
        }
        MyFile << endl;

        MyFile << "    " <<"CLIQUES FOUND: " << num_cliques_found << endl;

        




        now = chrono::system_clock::now();
       // Convert the current time to time since epoch
        duration2 = now.time_since_epoch();

        milliseconds2
            = chrono::duration_cast<chrono::nanoseconds>(
                  duration2)
                  .count();

        milliseconds = (double(milliseconds2)- double(milliseconds1))/1000000000;

       MyFile << "    " <<"TIME: " << milliseconds << endl;

       MyFile << "############################" << endl;

       iter++;

    }
}



void count_cliques(){

    std::ifstream infile("test_smiles");
    std::vector<std::string> smiles;
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        smiles.push_back(line);       
        
    }
    std::vector<std::string> smiles2;
    std::ifstream infile2("test_smiles");
    while (std::getline(infile2, line))
    {
        std::istringstream iss(line);
        smiles2.push_back(line);       
        
    }
    std::mt19937_64 random_engine{std::random_device{}()};

    int iter = 0;
    ofstream MyFile("clique_count.txt");
    while(iter<100000){
        std::vector<std::string> instance{};
        std::ranges::sample(smiles, std::back_inserter(instance), 4, random_engine);

        
        std::vector<LineGraph> graphs;
        auto g = smiles_to_graph(instance[0]);
        auto lg = molecular_graph_to_line_graph(g);
        graphs.push_back(lg);

        auto h = smiles_to_graph(instance[1]);
        auto hg = molecular_graph_to_line_graph(h);
        graphs.push_back(hg);
        
        auto q = smiles_to_graph(instance[2]);
        auto  qg = molecular_graph_to_line_graph(q);
        graphs.push_back(qg);
        
        

        std::vector<std::string> instance2{};
        std::ranges::sample(smiles2, std::back_inserter(instance2), 3, random_engine);

        auto w = smiles_to_graph(instance2[0]);
        auto  wg = molecular_graph_to_line_graph(w);
        graphs.push_back(wg);

        auto r = smiles_to_graph(instance2[1]);
        auto  rg = molecular_graph_to_line_graph(r);
        graphs.push_back(rg);

        auto rng = std::default_random_engine {};
        std::shuffle(std::begin(graphs), std::end(graphs), rng);

         auto now = chrono::system_clock::now();

        auto duration = now.time_since_epoch();

        auto milliseconds1
            = chrono::duration_cast<chrono::nanoseconds>(
                  duration)
                  .count();


       
        special_product = false;

        maximum_common_subgraph(graphs, true);


        MyFile << num_cliques_found << ",";


        now = chrono::system_clock::now();
       // Convert the current time to time since epoch
        auto duration2 = now.time_since_epoch();

        auto milliseconds2
            = chrono::duration_cast<chrono::nanoseconds>(
                  duration2)
                  .count();

        auto milliseconds = (double(milliseconds2)- double(milliseconds1))/1000000000;

       MyFile << milliseconds << endl;

       for(auto i : instance2){
            cout << i << endl;
       }
       for(auto i : instance){
            cout << i << endl;
       }

        cout  <<"MCS SIZE: " << maximum << endl;


       iter++;

    }



}




void test_time_ordering(){

    std::vector<std::vector<std::string>> instances;

     ofstream MyFile("mcs_results1.txt");

     ofstream results("test");

     //results << "solution_time,ordering_time" << endl;  

     

     double ordering_time = 0;
     

    std::ifstream infile("instances_wl3.txt"/*"random_smiles"*/);
     //std::ifstream infile("random_smiles");
    std::string line;
    std::vector<std::string> instance;
    while (std::getline(infile, line))
    {

        std::istringstream iss(line);

        if(line=="#"){
            
            auto v = instance;
            instances.push_back(v);
            instance = std::vector<std::string>();
            continue;
        }

        instance.push_back(line); 
          
        
    }

    cout << instances.size() << endl;
    

    

    int iter = 0;
    auto rng = std::default_random_engine {};
    rng.seed(50);

    
    for(auto &inst : instances){
        bool bad = false;
        auto instance = inst;
        std::string good_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCC";

        //cout << "[";
        for(auto i : instance){
            //cout << '\"' << i << '\"'<< ',' <<  endl;
            //cout << i << endl;    
            
            for(int j=0; j<i.size() ; j++){
                auto c = i[j];
                if(c=='+'){
                    bad = true;
                }
            }
            if(not bad){
                good_smiles = i;
            }

        }

        for(int k = 0; k<instance.size();k++){
            auto i = instance[k];
            //cout << '\"' << i << '\"'<< ',' <<  endl;
            //cout << i << endl; 
            bad = false;   
            
            for(int j=0; j<i.size() ; j++){
                auto c = i[j];
                if(c=='+'){
                   
                    bad = true;
                }
            }
            bad = false;
            if(bad){
                instance[k]=good_smiles;
            }

        }


        for(auto i : instance){
            //cout << '\"' << i << '\"'<< ',' <<  endl;
            //cout << i << endl;    
            
            for(int j=0; j<i.size() ; j++){
                auto c = i[j];
                if(c=='+'){
                  // cout << i << endl;
                }
            }

        }
        



        //cout << "]" << endl;

        


        
       
        
        /*for(auto i : instance){
            cout << "\"" <<  i << "\"," <<  endl;
        }*/
        
        /*for(auto k : instance){
            cout << k << endl;
        }*/
        if(iter==15000999){
            iter++;
            continue;
        }
        
        if(iter==2000){
            //break;
        }
        

              
        
        std::vector<LineGraph> graphs;
        std::vector<MolecularGraph> mgraphs;

        
        for(auto smiles : instance){

            
           

            auto g = smiles_to_graph(smiles);
            
             
            auto lgg =  molecular_graph_to_line_graph(g);

            graphs.push_back(lgg);
            

        }
        
        


            
        
        
        



        
       
        




    



        


       // std::shuffle(std::begin(graphs), std::end(graphs), rng);
        if(iter<=-1){//536//7/100/77
            iter++;
            continue;

        }

        for(auto j : instance){
            cout << "\"" << j << "\"," <<  endl;
        }


        auto now2 = chrono::system_clock::now();

            auto duration2 = now2.time_since_epoch();

            auto milliseconds3
                = chrono::duration_cast<chrono::nanoseconds>(
                      duration2)
                      .count();


        
        
        
        
        
         now2 = chrono::system_clock::now();
       // Convert the current time to time since epoch
        auto duration3 = now2.time_since_epoch();

        auto milliseconds4
            = chrono::duration_cast<chrono::nanoseconds>(
                  duration3)
                  .count();
        
        double time_to_order = (double(milliseconds4)- double(milliseconds3))/1000000000;
        
        
        ordering_time = ordering_time + time_to_order;

        
        //graphs = {graphs[0], graphs[1], graphs[2]};
       //std::shuffle(std::begin(graphs), std::end(graphs), rng);



        


        int num_nodes = 0;

        for(auto curr_graph : graphs){
            num_nodes = num_nodes + num_vertices(curr_graph.g);
        }

        double aver_size = num_nodes/graphs.size();

       
        special_product = false;

        auto now = chrono::system_clock::now();

        auto duration = now.time_since_epoch();

        auto milliseconds1
            = chrono::duration_cast<chrono::nanoseconds>(
                  duration)
                  .count();


        
        maximum_common_subgraph(graphs, true);


         now = chrono::system_clock::now();
       // Convert the current time to time since epoch
        auto duration5 = now.time_since_epoch();

        auto milliseconds2
            = chrono::duration_cast<chrono::nanoseconds>(
                  duration5)
                  .count();

        auto milliseconds = (double(milliseconds2)- double(milliseconds1))/1000000000;

       

        

        //cout << "################################" << endl;

        //cout << iter << ", " << num_edges(maximum_graph.g) << ", " << all_maximum_graphs.size() << ", " << num_cliques_found << endl;
        cout << "#########################################" << endl;
        cout << "INSTANCE: " << iter << endl;
        cout << "#GRAPHS: " << graphs.size() << endl;
        cout << "AVERAGE #EDGES: " << aver_size << endl;
        cout << "MCS SIZE: " << num_edges(maximum_graph.g) << endl;
        cout << "SOLUTIONS: " << endl;
        cout << "ORDERING TIME: " << time_to_order << endl;
        cout << "TIME: " << milliseconds << endl;
        for(auto mcs_graph : all_maximum_graphs){
            cout << to_mol(mcs_graph) << endl;
        }
        cout << "#########################################" << endl;
        cout << endl << endl;

        results << milliseconds + time_to_order << endl;

       

       
        int u = 0;

        MyFile << graphs.size() << ",";

        for(int j = 0; j<graphs.size(); j++){

            auto i = instance[j];

           
            MyFile << i << ",";
        }   
        u++;
        for(auto &t : all_maximum_graphs){
            std::string smi = to_mol(t);
            MyFile << smi << ",";
        }


        MyFile << endl;

        
    
        iter++;
        
        

    }

    cout << ordering_time << endl;

   
       




}

int clique_count_int = 0;





void generate_clique_finding_data(){


    std::vector<std::string> smiles;

    ofstream MyFile("clique_data.txt");
    std::ifstream infile("test_smiles");


    std::string line;

    while (std::getline(infile, line))
    {
        bool  bad = false;
        std::istringstream iss(line);
        for(auto i : line){
            if(i=='i'){
                bad = true;
            }
        }
        if(bad){
            continue;
        }
        smiles.push_back(line);
        cout << line << endl;       
        
    }




    int iter = 0;

    std::mt19937 generator(std::random_device{}());

    std::uniform_int_distribution<std::size_t> distribution(0, smiles.size() - 1);

    int prev = 0;

    special_product = true;

    while(iter<100000){
        iter ++;

        std::size_t index1 = distribution(generator);
        std::size_t index2 = distribution(generator);

        auto g1 = smiles_to_graph(smiles[index1]);
        auto lg1 = molecular_graph_to_line_graph(g1);
        
        auto g2 = smiles_to_graph(smiles[index2]);
        auto lg2 = molecular_graph_to_line_graph(g2);


        auto direct_product = modular_graph_product(lg1, lg2);

        auto type1_connected_components = find_type1_connected_components(direct_product);

        cout << "DONE!" << endl;

        auto now = chrono::system_clock::now();

        auto duration = now.time_since_epoch();

        auto milliseconds1
            = chrono::duration_cast<chrono::nanoseconds>(
                  duration)
                  .count();

        int cliques_num = 0;

        for(auto subgraph : type1_connected_components){


            

            auto cliques =  bron_kerbosch_driver(subgraph);
            //for(const auto& clique: cliques){
            while (cliques.move_next()){

                auto clique = cliques.current_value();

                cliques_num++;
            }

        }




        now = chrono::system_clock::now();
       // Convert the current time to time since epoch
        auto duration2 = now.time_since_epoch();

        auto milliseconds2
            = chrono::duration_cast<chrono::nanoseconds>(
                  duration2)
                  .count();

        auto milliseconds = (double(milliseconds2)- double(milliseconds1))/1000000000;

        cout << milliseconds << "," << cliques_num << "," <<smiles[index1] << "," << smiles[index2] << endl;
            
        MyFile << milliseconds << "," << cliques_num << "," <<smiles[index1] << "," << smiles[index2] << endl;






    }






}



void ordering_experiment(){

    std::vector<std::vector<std::string>> instances;

     ofstream MyFile("ordering_experiment_results.txt");

     cout << "LOL";
     

    std::ifstream infile("ordering_experiment_smiles.txt");
     //std::ifstream infile("random_smiles");
    std::string line;
    std::vector<std::string> instance;
    while (std::getline(infile, line))
    {

        std::istringstream iss(line);

        if(line=="#"){
            
            auto v = instance;
            instances.push_back(v);
            instance = std::vector<std::string>();
            continue;
        }

        instance.push_back(line);       
        
    }

    std::vector<std::vector<int>> permutations;

    std::vector<int> v = {0, 1,2,3, 4};

    do {
        //cout << v[0] << v[1] << v[2] << v[3] << v[4]<< endl;
        bool found = false;
        for(auto perm : permutations){
            if(perm[2]==v[2] and perm[3]==v[3] and perm[4]==v[4]){
                found = true;
            }
        }
        if(not found){
            permutations.push_back(v);
        }
        
    } while (std::next_permutation(v.begin(), v.end()));



    for(auto instance : instances){

        MyFile << instance[0] << endl << instance[1] << endl << instance[2] << endl << instance[3] << endl << instance[4] << endl;


        for(auto permutation : permutations){


            std::vector<LineGraph> graphs;

            auto g = smiles_to_graph(instance[permutation[0]]);
            auto lg = molecular_graph_to_line_graph(g);
            graphs.push_back(lg);

            auto h = smiles_to_graph(instance[permutation[1]]);
            auto hg = molecular_graph_to_line_graph(h);
            graphs.push_back(hg);

            
            
            auto q = smiles_to_graph(instance[permutation[2]]);
            auto  qg = molecular_graph_to_line_graph(q);
            graphs.push_back(qg);
            
            
            auto f = smiles_to_graph(instance[permutation[3]]);
            auto fg = molecular_graph_to_line_graph(f);
            graphs.push_back(fg);
            
            
            
            auto a = smiles_to_graph(instance[permutation[4]]);
            auto  ag = molecular_graph_to_line_graph(a);
            graphs.push_back(ag);



            auto now = chrono::system_clock::now();

            auto duration = now.time_since_epoch();

            auto milliseconds1
                = chrono::duration_cast<chrono::nanoseconds>(
                      duration)
                      .count();



            special_product = false;

            maximum_common_subgraph(graphs, true);


            MyFile << num_cliques_found << ",";





             now = chrono::system_clock::now();
           // Convert the current time to time since epoch
            auto duration2 = now.time_since_epoch();

            auto milliseconds2
                = chrono::duration_cast<chrono::nanoseconds>(
                      duration2)
                      .count();

            auto milliseconds = (double(milliseconds2)- double(milliseconds1))/1000000000;

            cout << "TIME: " << milliseconds << ", SIZE: " << num_edges(maximum_graph.g) << endl;


            MyFile << milliseconds << ",";

           MyFile << permutation[0] << permutation[1] << permutation[2] << permutation[3] << permutation[4] << endl;





        }
    }








}


/*******************************functions for experiments**********************************************************************/


int main(int argc, char* argv[]) {

    
    std::string input = argv[1];
    std::string output = argv[2];
    std::string product_type = argv[3];
    std::string ordering = argv[4];
    
    bool use_special_product = false;
    bool order = false;



    if(product_type=="prune"){
        use_special_product = true;
    }
    else if(product_type=="no-prune"){
        use_special_product = false;
    }
    else{
        cout << "Command line argument 3 should be \"prune\" or \"no-prune\"" << endl;
        return 0;
    }

    if(ordering=="order"){
        order = true;
    }
    else if(ordering=="no-order"){
        order = false;
    }
    else{
        cout << "Command line argument 4 should be \"order\" or \"no-order\"" << endl;
        return 0;
    }

    
    
    solve_instances(input, output, use_special_product, order);
    
    
    return 0;
}
