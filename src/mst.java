import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.TimeUnit;

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

class Node implements Comparable<Node> {
    int index, rank;
    double key;
    Node parent;

    // For Kruskal's Disjoint set
    public Node (int index, int rank, Node p){
        this.index = index;
        this.rank = rank;
        this.parent = p;
    }

    // For Prim's algorithm
    public Node (int index, double key, Node p){
        this.index = index;
        this.key = key;
        this.parent = p;
    }

    @Override
    public int compareTo(Node o) {
        if (this.key > o.key) return 1;
        else if (this.key < o.key) return -1;
        return 0;
    }
}

// Disjoint set implemented using LinkedList approach for Kruskal's algorithm
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
        Node current = n, next;
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
    ArrayList<Node> nodes;
    Boolean[] inQueue;

    public Graph (int n) {
        int m = n * (n - 1) / 2;    // number of edges in a complete graph
        adjList = new HashMap<Integer, List<Integer>>();
        vertices = new ArrayList<Vertex>(n);
        edges = new ArrayList<Edge>(m);
        inQueue = new Boolean[n];
    }

    /**
     *
     * @param n
     * @param complete
     * @return a connected or complete graph with n vertices
     */
    public static Graph randomGraph(int n, boolean complete) {
        Graph g = new Graph(n);

        // randomly generate vertices between (0,0) and (1,1)
        for (int i = 0; i < n; i++) {
            g.vertices.add(new Vertex(ThreadLocalRandom.current().nextDouble(),
                    ThreadLocalRandom.current().nextDouble()));
            g.adjList.put(i, new LinkedList<Integer>());
        }

        // randomly generate edges based on p
        Vertex u, v;
        double p;
        if (complete) p = 1.0;
        else p = 4 * Math.log(n) / n;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (ThreadLocalRandom.current().nextDouble() < p){
                    u = g.vertices.get(i);
                    v = g.vertices.get(j);
                    g.edges.add(new Edge(u, v));
                    g.adjList.get(i).add (j);
                    g.adjList.get(j).add (i);
                }
            }
        }
        return g;
    }

    public void createDisjointSet(){
        ds = new DisjointSet(vertices);
    }

    public void initializePrim(){
        nodes = new ArrayList<Node>(vertices.size());
        for (int i = 0; i < vertices.size(); i++){
            Node n = new Node (nodes.size(), Double.MAX_VALUE, null);
            inQueue[i] = true;
            vertices.get(i).n = n;
            nodes.add(n);
        }
        inQueue[0] = false;
        nodes.get(0).key = 0;

    }
}

public class mst {

    static void runKruskalComplete (int n, int p){
        ArrayList<Edge> mst;
        Graph g;
        double weight = 0;
        for (int i = 0; i < p; i++){
            mst = new ArrayList<>();
            g = Graph.randomGraph(n, true);
            g.createDisjointSet();
            Collections.sort(g.edges);
            for (Edge e : g.edges){
                Vertex u = e.u;
                Vertex v = e.v;
                if (g.ds.findSet(u.n) != g.ds.findSet(v.n)){
                    weight += e.w;
                    mst.add (e);
                    g.ds.union(u.n, v.n);
                    // MST found
                    if (g.ds.sets == 1){
                        break;
                    }
                }
            }
        }
        System.out.printf("%6f\n", weight/50.0);
    }

    static void runKruskal (int n, int p){
        ArrayList<Edge> mst;
        Graph g;
        long runTime = 0;
        for (int i = 0; i < p; i++){
            mst = new ArrayList<>();
            g = Graph.randomGraph(n, false);
            long startTime = System.currentTimeMillis();
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
            runTime += (System.currentTimeMillis() - startTime);
        }
        System.out.print((runTime / 50.0) +"ms and ");
    }

    static void runPrim (int n, int p){
        long runTime = 0;
        Graph g;
        PriorityQueue<Node> queue;
        for (int pi = 0; pi < p; pi++){
            g = Graph.randomGraph(n, false);
            long startTime = System.currentTimeMillis();
            g.initializePrim();
            queue = new PriorityQueue<>(n);

            for (int i = 0; i < n; i++) {
                queue.add(g.nodes.get(i));
            }
            int count = 0;
            while (!queue.isEmpty()){
                Node current = queue.poll();
                Vertex u = g.vertices.get (current.index);
                for (Integer d : g.adjList.get(current.index)){
                    if (g.inQueue[d]){
                        Vertex v = g.vertices.get (d);
                        if (Vertex.dist(u,v) < v.n.key) {
                            queue.remove(v.n);
                            v.n.parent = u.n;
                            v.n.key = Vertex.dist(u, v);
                            queue.add(v.n);
                        }
                    }
                }
            }
            runTime += System.currentTimeMillis() - startTime;
        }
        System.out.println((runTime / 50.0) +"ms respectively.");
    }


    public static void main(String[] args) throws InterruptedException {

        System.out.print("When n = 100:  ");
        System.out.print("The average weight of MSTs in 50 complete graphs is ");
        runKruskalComplete(100,50);

        System.out.print("When n = 500:  ");
        System.out.print("The average weight of MSTs in 50 complete graphs is ");
        runKruskalComplete(500,50);

        System.out.print("When n = 1000: ");
        System.out.print("The average weight of MSTs in 50 complete graphs is ");
        runKruskalComplete(1000,50);

        System.out.print("When n = 5000: ");
        System.out.print("The average weight of MSTs in 50 complete graphs is ");
        runKruskalComplete(5000,50);

        TimeUnit.SECONDS.sleep(1);

        System.out.print("When n = 100:  The average run time for Kruskal and Prim on 50 connected graphs are " );
        runKruskal (100, 50);
        runPrim (100, 50);

        System.out.print("When n = 500:  The average run time for Kruskal and Prim on 50 connected graphs are " );
        runKruskal (500, 50);
        runPrim (500, 50);

        System.out.print("When n = 1000: The average run time for Kruskal and Prim on 50 connected graphs are " );
        runKruskal (1000, 50);
        runPrim (1000, 50);

        System.out.print("When n = 5000: The average run time for Kruskal and Prim on 50 connected graphs are " );
        runKruskal (5000, 50);
        runPrim (5000, 50);
    }
}
