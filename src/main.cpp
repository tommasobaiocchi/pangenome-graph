#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <stdexcept>
#include <functional>

using namespace std;

class Graph {
public:
    struct Node {
        int id;
        string sequence;
        vector<int> forward_edges;   // archi in direzione +
        vector<int> reverse_edges;   // archi in direzione -
        
        Node(int node_id, const string& seq) : id(node_id), sequence(seq) {}
    };

private:
    unordered_map<int, shared_ptr<Node>> nodes;
    
    bool dfsCycleDetection(int node_id, vector<bool>& visited, vector<bool>& in_stack, 
                          vector<pair<int, int>>& back_edges) const {
        int original_id = node_id / 2;
        bool is_reverse = (node_id % 2 == 1);
        
        if (!visited[node_id]) {
            visited[node_id] = true;
            in_stack[node_id] = true;
            
            const auto& neighbors = is_reverse ? 
                nodes.at(original_id)->reverse_edges : 
                nodes.at(original_id)->forward_edges;
                
            for (int neighbor : neighbors) {
                int neighbor_id = (neighbor > 0) ? neighbor * 2 : (-neighbor * 2 + 1);
                
                if (!visited[neighbor_id]) {
                    if (dfsCycleDetection(neighbor_id, visited, in_stack, back_edges)) {
                        return true;
                    }
                } else if (in_stack[neighbor_id]) {
                    back_edges.push_back({node_id, neighbor_id});
                }
            }
        }
        
        in_stack[node_id] = false;
        return false;
    }

public:
    Graph() {}
    
    void addNode(int id, const string& sequence) {
        if (nodes.find(id) == nodes.end()) {
            nodes[id] = make_shared<Node>(id, sequence);
        } else {
            nodes[id]->sequence = sequence;
        }
    }
    
    void addEdge(int from, char from_orient, int to, char to_orient) {
        if (nodes.find(from) == nodes.end() || nodes.find(to) == nodes.end()) {
            throw invalid_argument("Node ID not found in graph");
        }
        
        int target_node = to;
        if (to_orient == '-') {
            target_node = -to;
        }
        
        if (from_orient == '+') {
            nodes[from]->forward_edges.push_back(target_node);
        } else {
            nodes[from]->reverse_edges.push_back(target_node);
        }
    }
    
    const Node* getNode(int id) const {
        auto it = nodes.find(id);
        return (it != nodes.end()) ? it->second.get() : nullptr;
    }
    
    vector<pair<int, int>> findBackEdges() const {
        vector<pair<int, int>> back_edges;
        if (nodes.empty()) return back_edges;
        
        vector<bool> visited(nodes.size() * 2, false);
        vector<bool> in_stack(nodes.size() * 2, false);
        
        for (const auto& pair : nodes) {
            int node_id = pair.first;
            if (!visited[node_id * 2]) {
                dfsCycleDetection(node_id * 2, visited, in_stack, back_edges);
            }
        }
        
        return back_edges;
    }
    
    void removeBackEdges() {
        auto back_edges = findBackEdges();
        for (const auto& edge : back_edges) {
            int from_id = edge.first / 2;
            bool from_reverse = (edge.first % 2 == 1);
            int to_id = edge.second / 2;
            
            auto& edges = from_reverse ? nodes[from_id]->reverse_edges : nodes[from_id]->forward_edges;
            int target = to_id * 2;
            if (edge.second % 2 == 1) target = -to_id;
            
            edges.erase(remove(edges.begin(), edges.end(), target), edges.end());
        }
    }
    
    bool isAcyclic() const {
        return findBackEdges().empty();
    }
    
    vector<int> findSources() const {
        unordered_map<int, int> in_degree;
        
        for (const auto& pair : nodes) {
            in_degree[pair.first] = 0;
        }
        
        for (const auto& pair : nodes) {
            for (int neighbor : pair.second->forward_edges) {
                int target_id = (neighbor > 0) ? neighbor : -neighbor;
                in_degree[target_id]++;
            }
            for (int neighbor : pair.second->reverse_edges) {
                int target_id = (neighbor > 0) ? neighbor : -neighbor;
                in_degree[target_id]++;
            }
        }
        
        vector<int> sources;
        for (const auto& pair : in_degree) {
            if (pair.second == 0) {
                sources.push_back(pair.first);
            }
        }
        
        return sources;
    }
    
    vector<int> findSinks() const {
        vector<int> sinks;
        
        for (const auto& pair : nodes) {
            const auto& node = pair.second;
            if (node->forward_edges.empty() && node->reverse_edges.empty()) {
                sinks.push_back(node->id);
            }
        }
        
        return sinks;
    }
    
    size_t getNodeCount() const { return nodes.size(); }
    const unordered_map<int, shared_ptr<Node>>& getNodes() const { return nodes; }
};

class RollingHash {
private:
    const long long base = 911382323;
    const long long prime = 1000000007;
    long long base_power;
    int k;
    
    long long charToValue(char c) const {
        switch(c) {
            case 'A': return 1;
            case 'T': return 2;
            case 'C': return 3;
            case 'G': return 4;
            default: return 0;
        }
    }

public:
    RollingHash(int kmer_length) : k(kmer_length) {
        base_power = 1;
        for (int i = 0; i < k-1; i++) {
            base_power = (base_power * base) % prime;
        }
    }
    
    long long compute(const string& s) const {
        if (s.length() != k) {
            throw invalid_argument("String length must equal k");
        }
        
        long long hash = 0;
        for (char c : s) {
            hash = (hash * base + charToValue(c)) % prime;
        }
        return hash;
    }
    
    long long update(long long old_hash, char remove_char, char add_char) const {
        long long remove_val = charToValue(remove_char);
        long long add_val = charToValue(add_char);
        
        long long hash = (old_hash - remove_val * base_power) % prime;
        hash = (hash * base + add_val) % prime;
        
        if (hash < 0) hash += prime;
        return hash;
    }
};

class PathFinder {
public:
    struct OrientedNode {
        int id;
        bool is_reverse;
        
        OrientedNode(int node_id, bool reverse = false) : id(node_id), is_reverse(reverse) {}
        
        bool operator==(const OrientedNode& other) const {
            return id == other.id && is_reverse == other.is_reverse;
        }
    };
    
    struct Path {
        vector<OrientedNode> nodes;
        string sequence;
    };
    
private:
    const Graph& graph;
    
    string reverseComplement(const string& sequence) const {
        string result;
        result.reserve(sequence.length());
        
        for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
            switch (*it) {
                case 'A': result += 'T'; break;
                case 'T': result += 'A'; break;
                case 'C': result += 'G'; break;
                case 'G': result += 'C'; break;
                default: result += *it; break;
            }
        }
        
        return result;
    }
    
    string getNodeSequence(const OrientedNode& node) const {
        const auto* graph_node = graph.getNode(node.id);
        if (!graph_node) return "";
        
        return node.is_reverse ? reverseComplement(graph_node->sequence) : graph_node->sequence;
    }
    
    void dfsFindPaths(const OrientedNode& current, const OrientedNode& target, 
                     vector<OrientedNode>& current_path, vector<Path>& all_paths,
                     int max_paths = 1000) const {
        if (all_paths.size() >= max_paths) return;
        
        current_path.push_back(current);
        
        if (current.id == target.id && current.is_reverse == target.is_reverse) {
            Path path;
            path.nodes = current_path;
            path.sequence = buildSequenceFromPath(current_path);
            all_paths.push_back(path);
        } else {
            const auto* graph_node = graph.getNode(current.id);
            if (graph_node) {
                // Esplora archi forward
                for (int neighbor : graph_node->forward_edges) {
                    OrientedNode next(neighbor > 0 ? neighbor : -neighbor, neighbor < 0);
                    if (find(current_path.begin(), current_path.end(), next) == current_path.end()) {
                        dfsFindPaths(next, target, current_path, all_paths, max_paths);
                    }
                }
                
                // Esplora archi reverse
                for (int neighbor : graph_node->reverse_edges) {
                    OrientedNode next(neighbor > 0 ? neighbor : -neighbor, neighbor < 0);
                    if (find(current_path.begin(), current_path.end(), next) == current_path.end()) {
                        dfsFindPaths(next, target, current_path, all_paths, max_paths);
                    }
                }
            }
        }
        
        current_path.pop_back();
    }
    
    string buildSequenceFromPath(const vector<OrientedNode>& path) const {
        string sequence;
        
        for (const auto& node : path) {
            sequence += getNodeSequence(node);
        }
        
        return sequence;
    }

public:
    PathFinder(const Graph& graph_ref) : graph(graph_ref) {}
    
    vector<Path> findAllPaths(int source, int target, int max_paths = 1000) const {
        vector<Path> all_paths;
        vector<OrientedNode> current_path;
        
        dfsFindPaths(OrientedNode(source, false), OrientedNode(target, false), 
                    current_path, all_paths, max_paths);
        return all_paths;
    }
    
    vector<string> getAllPathSequences(int source, int target, int max_paths = 1000) const {
        auto paths = findAllPaths(source, target, max_paths);
        vector<string> sequences;
        
        for (const auto& path : paths) {
            sequences.push_back(path.sequence);
        }
        
        return sequences;
    }
};

class KmerAnalyzer {
private:
    const Graph& graph;
    int k;
    RollingHash hasher;
    PathFinder path_finder;
    
public:
    KmerAnalyzer(const Graph& graph_ref, int kmer_length) 
        : graph(graph_ref), k(kmer_length), hasher(kmer_length), path_finder(graph_ref) {}
    
    unordered_map<string, int> countKmersWithHash(const string& sequence) const {
        unordered_map<string, int> kmer_counts;
        if (sequence.length() < k) return kmer_counts;
        
        string first_kmer = sequence.substr(0, k);
        long long current_hash = hasher.compute(first_kmer);
        kmer_counts[first_kmer]++;
        
        for (size_t i = 1; i <= sequence.length() - k; i++) {
            current_hash = hasher.update(current_hash, sequence[i-1], sequence[i+k-1]);
            
            // Per sicurezza, verifica con confronto diretto (in produzione si userebbero multiple hash functions)
            string kmer = sequence.substr(i, k);
            kmer_counts[kmer]++;
        }
        
        return kmer_counts;
    }
    
    unordered_map<string, int> countKmers(const string& sequence) const {
        return countKmersWithHash(sequence);
    }
    
    vector<pair<string, int>> getTopKmers(const string& sequence, int top_n) const {
        auto counts = countKmers(sequence);
        vector<pair<string, int>> sorted_counts(counts.begin(), counts.end());
        
        sort(sorted_counts.begin(), sorted_counts.end(),
              [](const auto& a, const auto& b) {
                  return a.second > b.second || 
                        (a.second == b.second && a.first < b.first);
              });
        
        if (top_n > 0 && top_n < sorted_counts.size()) {
            sorted_counts.resize(top_n);
        }
        
        return sorted_counts;
    }
    
    unordered_map<string, int> countKmersInAllPathsEfficient(int source, int target, int max_paths = 100) const {
        unordered_map<string, int> total_counts;
        auto paths = path_finder.findAllPaths(source, target, max_paths);
        
        for (const auto& path : paths) {
            auto counts = countKmersWithHash(path.sequence);
            for (const auto& pair : counts) {
                total_counts[pair.first] += pair.second;
            }
        }
        
        return total_counts;
    }
    
    unordered_map<string, int> countKmersInAllPaths(int source, int target) const {
        return countKmersInAllPathsEfficient(source, target, 100);
    }
    
    vector<pair<string, int>> getTopKmersInAllPaths(int source, int target, int top_n) const {
        auto total_counts = countKmersInAllPaths(source, target);
        vector<pair<string, int>> sorted_counts(total_counts.begin(), total_counts.end());
        
        sort(sorted_counts.begin(), sorted_counts.end(),
              [](const auto& a, const auto& b) {
                  return a.second > b.second || 
                        (a.second == b.second && a.first < b.first);
              });
        
        if (top_n > 0 && top_n < sorted_counts.size()) {
            sorted_counts.resize(top_n);
        }
        
        return sorted_counts;
    }
    
    // Pattern matching con rolling hash
    vector<int> findPatternWithHash(const string& sequence, const string& pattern) const {
        vector<int> positions;
        int m = pattern.length();
        if (m == 0 || sequence.length() < m) return positions;
        
        if (m == k) {
            // Usa rolling hash per pattern della stessa lunghezza k
            long long pattern_hash = hasher.compute(pattern);
            string first_kmer = sequence.substr(0, m);
            long long current_hash = hasher.compute(first_kmer);
            
            if (current_hash == pattern_hash && first_kmer == pattern) {
                positions.push_back(0);
            }
            
            for (size_t i = 1; i <= sequence.length() - m; i++) {
                current_hash = hasher.update(current_hash, sequence[i-1], sequence[i+m-1]);
                if (current_hash == pattern_hash) {
                    // Verifica per evitare collisioni
                    if (sequence.substr(i, m) == pattern) {
                        positions.push_back(i);
                    }
                }
            }
        } else {
            // Per pattern di lunghezza diversa, usa ricerca diretta
            for (size_t i = 0; i <= sequence.length() - m; i++) {
                if (sequence.compare(i, m, pattern) == 0) {
                    positions.push_back(i);
                }
            }
        }
        
        return positions;
    }
    
    vector<int> findPattern(const string& sequence, const string& pattern) const {
        return findPatternWithHash(sequence, pattern);
    }
    
    vector<pair<int, int>> findPatternInAllPaths(int source, int target, const string& pattern) const {
        vector<pair<int, int>> results;
        auto paths = path_finder.findAllPaths(source, target, 100);
        
        for (size_t path_idx = 0; path_idx < paths.size(); path_idx++) {
            auto positions = findPatternWithHash(paths[path_idx].sequence, pattern);
            for (int pos : positions) {
                results.push_back({static_cast<int>(path_idx), pos});
            }
        }
        
        return results;
    }
};

class GFAParser {
private:
    static void parseSegment(const string& line, Graph& graph) {
        stringstream ss(line);
        char record_type;
        int node_id;
        string sequence_field;
        
        ss >> record_type >> node_id >> sequence_field;
        
        if (ss.fail()) {
            throw invalid_argument("Malformed segment line");
        }
        
        graph.addNode(node_id, sequence_field);
    }
    
    static void parseLink(const string& line, Graph& graph) {
        stringstream ss(line);
        char record_type;
        int from_id, to_id;
        char from_orient, to_orient;
        string overlap;
        
        ss >> record_type >> from_id >> from_orient >> to_id >> to_orient >> overlap;
        
        if (ss.fail() || (from_orient != '+' && from_orient != '-') || 
            (to_orient != '+' && to_orient != '-')) {
            throw invalid_argument("Malformed link line");
        }
        
        graph.addEdge(from_id, from_orient, to_id, to_orient);
    }
    
    static void processLine(const string& line, Graph& graph) {
        stringstream ss(line);
        char record_type;
        ss >> record_type;
        
        switch (record_type) {
            case 'S':
                parseSegment(line, graph);
                break;
            case 'L':
                parseLink(line, graph);
                break;
            default:
                break;
        }
    }

public:
    static unique_ptr<Graph> parseFile(const string& filename) {
        auto graph = make_unique<Graph>();
        
        ifstream file(filename);
        if (!file.is_open()) {
            throw runtime_error("Cannot open GFA file: " + filename);
        }
        
        string line;
        int line_number = 0;
        
        while (getline(file, line)) {
            line_number++;
            if (line.empty() || line[0] == '#') continue;
            
            try {
                processLine(line, *graph);
            } catch (const exception& e) {
                cerr << "Warning: Error parsing line " << line_number << ": " << e.what() << endl;
            }
        }
        
        return graph;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <gfa_file> [kmer_length] [top_n] [source] [target]" << endl;
        cout << "Example: " << argv[0] << " example.gfa 9 10 1 100" << endl;
        return 1;
    }
    
    string filename = argv[1];
    int kmer_length = (argc > 2) ? stoi(argv[2]) : 9;
    int top_n = (argc > 3) ? stoi(argv[3]) : 10;
    int source_id = (argc > 4) ? stoi(argv[4]) : -1;
    int target_id = (argc > 5) ? stoi(argv[5]) : -1;
    
    try {
        cout << "Parsing GFA file: " << filename << endl;
        auto graph = GFAParser::parseFile(filename);
        cout << "Loaded graph with " << graph->getNodeCount() << " nodes" << endl;
        
        // Rimuovi back edges per ottenere DAG
        auto back_edges = graph->findBackEdges();
        if (!back_edges.empty()) {
            cout << "Found " << back_edges.size() << " back edges. Removing them..." << endl;
            graph->removeBackEdges();
            cout << "Graph is now acyclic: " << (graph->isAcyclic() ? "Yes" : "No") << endl;
        } else {
            cout << "Graph is already acyclic" << endl;
        }
        
        auto sources = graph->findSources();
        auto sinks = graph->findSinks();
        
        cout << "Found " << sources.size() << " sources and " << sinks.size() << " sinks" << endl;
        
        if (sources.empty() || sinks.empty()) {
            cout << "Error: No sources or sinks found" << endl;
            return 1;
        }
        
        if (source_id == -1) source_id = sources[0];
        if (target_id == -1) target_id = sinks[0];
        
        cout << "Using source: " << source_id << ", target: " << target_id << endl;
        
        if (!graph->getNode(source_id)) {
            cout << "Error: Source node " << source_id << " not found" << endl;
            return 1;
        }
        if (!graph->getNode(target_id)) {
            cout << "Error: Target node " << target_id << " not found" << endl;
            return 1;
        }
        
        KmerAnalyzer analyzer(*graph, kmer_length);
        
        cout << "Finding top " << top_n << " " << kmer_length << "-mers..." << endl;
        auto top_kmers = analyzer.getTopKmersInAllPaths(source_id, target_id, top_n);
        
        cout << "Top " << top_n << " " << kmer_length << "-mers:" << endl;
        for (size_t i = 0; i < top_kmers.size(); i++) {
            cout << i + 1 << ". " << top_kmers[i].first << " : " << top_kmers[i].second << " occurrences" << endl;
        }
        
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}