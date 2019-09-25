import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

class Vertex {
    double x;
    double y;
    Node n;

    public Vertex(double x, double y) {
        this.x = x;
        this.y = y;
    }

    // Returns the Euclidean distance between two vertices
    public static double dist (Vertex a, Vertex b){
        return Math.sqrt((b.y-a.y)*(b.y-a.y) + (b.x-a.x)*(b.x-a.x));
    }

}

class Edge implements Comparable<Edge> {
    Vertex u, v;
    double w;

    public Edge(Vertex u, Vertex v) {
        this.u = u;
        this.v = v;
        this.w = Vertex.dist (u,v);
    }

    @Override
    public int compareTo(Edge o) {
        return Double.compare (this.w, o.w);
    }
}

// Node in a Disjoint Set
class Node {
    int index, rank;
    Node parent;

    public Node (int index, int rank, Node p){
        this.index = index;
        this.rank = rank;
        parent = p;
    }
}

// Disjoint set implemented using LinkedList approach
class DisjointSet{
    int sets = 0;
    ArrayList<Node> rootNodes;

    public void makeSet(Vertex vertex){
        Node n = new Node (rootNodes.size(), 0, null);
        vertex.n = n;
        rootNodes.add (n);
        sets++;
    }

    public int findSet(Node n){
        Node root = n;
        while (root.parent != null){
            root = root.parent;
        }

        // Path compression
        Node current = n, next = n.parent;
        while (current != root){
            next = current.parent;
            current.parent = root;
            current = next;
        }

        return root.index;
    }

    public void union(Node a, Node b){
        int aSetRoot = findSet(a);
        int bSetRoot = findSet(b);
        if (aSetRoot == bSetRoot) return;
        a = rootNodes.get(aSetRoot);
        b = rootNodes.get(bSetRoot);
        // Union by rank
        if (a.rank < b.rank){
            a.parent = b;
        }
        else if (a.rank > b.rank){
            b.parent = a;
        }
        else{
            a.parent = b;
            b.rank++;
        }
        sets--;
    }

    public DisjointSet (List<Vertex> vertices){
        rootNodes = new ArrayList<Node>(vertices.size());
        for (Vertex v : vertices){
            makeSet(v);
        }
    }
}
class Graph{
    Map<Integer, List<Integer>> adjList;
    ArrayList<Vertex> vertices;
    ArrayList<Edge> edges;
    DisjointSet ds;

    public Graph (int v) {
        int edgeSize = v * (v - 1) / 2;    // number of edges in a complete graph
        adjList = new HashMap<Integer, List<Integer>>();
        vertices = new ArrayList<Vertex>(v);
        edges = new ArrayList<Edge>(edgeSize);
    }

    /**
     *
     * @param n : number of vertices
     * @return a random complete graph with n vertices
     */
    public static Graph randomCompleteGraph(int n) {
        Graph g = new Graph(n);

        // randomly generate vertices between (0,0) and (1,1)
        for (int i = 0; i < n; i++) {
            g.vertices.add(new Vertex(ThreadLocalRandom.current().nextDouble(),
                    ThreadLocalRandom.current().nextDouble()));
            g.adjList.put(i, new LinkedList<Integer>());
        }

        Vertex a, b;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                a = g.vertices.get(i);
                b = g.vertices.get(j);
                g.edges.add(new Edge(a, b));
                g.adjList.get(i).add (j);
                g.adjList.get(j).add (i);
            }
        }
        return g;
    }

    public void createDisjointSet(){
        ds = new DisjointSet(vertices);
    }
}
public class Kruskal_Complete {

    static void runKruskal (int n, int p){
        double sum = 0;
        ArrayList<Edge> mst;
        Graph g;

        for (int i = 0; i < p; i++){
            mst = new ArrayList<>();
            g = Graph.randomCompleteGraph(n);
            g.createDisjointSet();
            Collections.sort(g.edges);
            for (Edge e : g.edges){
                Vertex u = e.u;
                Vertex v = e.v;
                if (g.ds.findSet(u.n) != g.ds.findSet(v.n)){
                    mst.add (e);
                    g.ds.union(u.n, v.n);
                    // MST found
                    if (g.ds.sets == 1){
                        break;
                    }
                }
            }
            for (Edge e : mst){
                sum += e.w;
            }
        }
        System.out.println("Average weight of MST for " + p + " complete graphs with " + n + " vertices: " +sum/p);
    }


    public static void main(String[] args) {
        runKruskal (100, 50);
        runKruskal (500, 50);
        runKruskal (1000, 50);
        runKruskal (5000, 50);
    }
}
