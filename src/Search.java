import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Search {
    private static final long INF = 100000000000000L;
    private static List<Vertex> VERTICES =
            Arrays.asList(
                    new Vertex("albanyGA",        31.58,  84.17),
                    new Vertex("albanyNY",        42.66,  73.78),
                    new Vertex("albuquerque",     35.11, 106.61),
                    new Vertex("atlanta",         33.76,  84.40),
                    new Vertex("augusta",         33.43,  82.02),
                    new Vertex("austin",          30.30,  97.75),
                    new Vertex("bakersfield",     35.36, 119.03),
                    new Vertex("baltimore",       39.31,  76.62),
                    new Vertex("batonRouge",      30.46,  91.14),
                    new Vertex("beaumont",        30.08,  94.13),
                    new Vertex("boise",           43.61, 116.24),
                    new Vertex("boston",          42.32,  71.09),
                    new Vertex("buffalo",         42.90,  78.85),
                    new Vertex("calgary",         51.00, 114.00),
                    new Vertex("charlotte",       35.21,  80.83),
                    new Vertex("chattanooga",     35.05,  85.27),
                    new Vertex("chicago",         41.84,  87.68),
                    new Vertex("cincinnati",      39.14,  84.50),
                    new Vertex("cleveland",       41.48,  81.67),
                    new Vertex("coloradoSprings", 38.86, 104.79),
                    new Vertex("columbus",        39.99,  82.99),
                    new Vertex("dallas",          32.80,  96.79),
                    new Vertex("dayton",          39.76,  84.20),
                    new Vertex("daytonaBeach",    29.21,  81.04),
                    new Vertex("denver",          39.73, 104.97),
                    new Vertex("desMoines",       41.59,  93.62),
                    new Vertex("elPaso",          31.79, 106.42),
                    new Vertex("eugene",          44.06, 123.11),
                    new Vertex("europe",          48.87,  -2.33),
                    new Vertex("ftWorth",         32.74,  97.33),
                    new Vertex("fresno",          36.78, 119.79),
                    new Vertex("grandJunction",   39.08, 108.56),
                    new Vertex("greenBay",        44.51,  88.02),
                    new Vertex("greensboro",      36.08,  79.82),
                    new Vertex("houston",         29.76,  95.38),
                    new Vertex("indianapolis",    39.79,  86.15),
                    new Vertex("jacksonville",    30.32,  81.66),
                    new Vertex("japan",           35.68, 220.23),
                    new Vertex("kansasCity",      39.08,  94.56),
                    new Vertex("keyWest",         24.56,  81.78),
                    new Vertex("lafayette",       30.21,  92.03),
                    new Vertex("lakeCity",        30.19,  82.64),
                    new Vertex("laredo",          27.52,  99.49),
                    new Vertex("lasVegas",        36.19, 115.22),
                    new Vertex("lincoln",         40.81,  96.68),
                    new Vertex("littleRock",      34.74,  92.33),
                    new Vertex("losAngeles",      34.03, 118.17),
                    new Vertex("macon",           32.83,  83.65),
                    new Vertex("medford",         42.33, 122.86),
                    new Vertex("memphis",         35.12,  89.97),
                    new Vertex("mexia",           31.68,  96.48),
                    new Vertex("mexico",          19.40,  99.12),
                    new Vertex("miami",           25.79,  80.22),
                    new Vertex("midland",         43.62,  84.23),
                    new Vertex("milwaukee",       43.05,  87.96),
                    new Vertex("minneapolis",     44.96,  93.27),
                    new Vertex("modesto",         37.66, 120.99),
                    new Vertex("montreal",        45.50,  73.67),
                    new Vertex("nashville",       36.15,  86.76),
                    new Vertex("newHaven",        41.31,  72.92),
                    new Vertex("newOrleans",      29.97,  90.06),
                    new Vertex("newYork",         40.70,  73.92),
                    new Vertex("norfolk",         36.89,  76.26),
                    new Vertex("oakland",         37.80, 122.23),
                    new Vertex("oklahomaCity",    35.48,  97.53),
                    new Vertex("omaha",           41.26,  96.01),
                    new Vertex("orlando",         28.53,  81.38),
                    new Vertex("ottawa",          45.42,  75.69),
                    new Vertex("pensacola",       30.44,  87.21),
                    new Vertex("philadelphia",    40.72,  76.12),
                    new Vertex("phoenix",         33.53, 112.08),
                    new Vertex("pittsburgh",      40.40,  79.84),
                    new Vertex("pointReyes",      38.07, 122.81),
                    new Vertex("portland",        45.52, 122.64),
                    new Vertex("providence",      41.80,  71.36),
                    new Vertex("provo",           40.24, 111.66),
                    new Vertex("raleigh",         35.82,  78.64),
                    new Vertex("redding",         40.58, 122.37),
                    new Vertex("reno",            39.53, 119.82),
                    new Vertex("richmond",        37.54,  77.46),
                    new Vertex("rochester",       43.17,  77.61),
                    new Vertex("sacramento",      38.56, 121.47),
                    new Vertex("salem",           44.93, 123.03),
                    new Vertex("salinas",         36.68, 121.64),
                    new Vertex("saltLakeCity",    40.75, 111.89),
                    new Vertex("sanAntonio",      29.45,  98.51),
                    new Vertex("sanDiego",        32.78, 117.15),
                    new Vertex("sanFrancisco",    37.76, 122.44),
                    new Vertex("sanJose",         37.30, 121.87),
                    new Vertex("sanLuisObispo",   35.27, 120.66),
                    new Vertex("santaFe",         35.67, 105.96),
                    new Vertex("saultSteMarie",   46.49,  84.35),
                    new Vertex("savannah",        32.05,  81.10),
                    new Vertex("seattle",         47.63, 122.33),
                    new Vertex("stLouis",         38.63,  90.24),
                    new Vertex("stamford",        41.07,  73.54),
                    new Vertex("stockton",        37.98, 121.30),
                    new Vertex("tallahassee",     30.45,  84.27),
                    new Vertex("tampa",           27.97,  82.46),
                    new Vertex("thunderBay",      48.38,  89.25),
                    new Vertex("toledo",          41.67,  83.58),
                    new Vertex("toronto",         43.65,  79.38),
                    new Vertex("tucson",          32.21, 110.92),
                    new Vertex("tulsa",           36.13,  95.94),
                    new Vertex("uk1",             51.30,   0.00),
                    new Vertex("uk2",             51.30,   0.00),
                    new Vertex("vancouver",       49.25, 123.10),
                    new Vertex("washington",      38.91,  77.01),
                    new Vertex("westPalmBeach",   26.71,  80.05),
                    new Vertex("wichita",         37.69,  97.34),
                    new Vertex("winnipeg",        49.90,  97.13),
                    new Vertex("yuma",            32.69, 114.62)
            );

    private static List<EdgeData> EDGES =
            Arrays.asList(
                    new EdgeData("albanyNY", "montreal", 226),
                    new EdgeData("albanyNY", "boston", 166),
                    new EdgeData("albanyNY", "rochester", 148),
                    new EdgeData("albanyGA", "tallahassee", 120),
                    new EdgeData("albanyGA", "macon", 106),
                    new EdgeData("albuquerque", "elPaso", 267),
                    new EdgeData("albuquerque", "santaFe", 61),
                    new EdgeData("atlanta", "macon", 82),
                    new EdgeData("atlanta", "chattanooga", 117),
                    new EdgeData("augusta", "charlotte", 161),
                    new EdgeData("augusta", "savannah", 131),
                    new EdgeData("austin", "houston", 186),
                    new EdgeData("austin", "sanAntonio", 79),
                    new EdgeData("bakersfield", "losAngeles", 112),
                    new EdgeData("bakersfield", "fresno", 107),
                    new EdgeData("baltimore", "philadelphia", 102),
                    new EdgeData("baltimore", "washington", 45),
                    new EdgeData("batonRouge", "lafayette", 50),
                    new EdgeData("batonRouge", "newOrleans", 80),
                    new EdgeData("beaumont", "houston", 69),
                    new EdgeData("beaumont", "lafayette", 122),
                    new EdgeData("boise", "saltLakeCity", 349),
                    new EdgeData("boise", "portland", 428),
                    new EdgeData("boston", "providence", 51),
                    new EdgeData("buffalo", "toronto", 105),
                    new EdgeData("buffalo", "rochester", 64),
                    new EdgeData("buffalo", "cleveland", 191),
                    new EdgeData("calgary", "vancouver", 605),
                    new EdgeData("calgary", "winnipeg", 829),
                    new EdgeData("charlotte", "greensboro", 91),
                    new EdgeData("chattanooga", "nashville", 129),
                    new EdgeData("chicago", "milwaukee", 90),
                    new EdgeData("chicago", "midland", 279),
                    new EdgeData("cincinnati", "indianapolis", 110),
                    new EdgeData("cincinnati", "dayton", 56),
                    new EdgeData("cleveland", "pittsburgh", 157),
                    new EdgeData("cleveland", "columbus", 142),
                    new EdgeData("coloradoSprings", "denver", 70),
                    new EdgeData("coloradoSprings", "santaFe", 316),
                    new EdgeData("columbus", "dayton", 72),
                    new EdgeData("dallas", "denver", 792),
                    new EdgeData("dallas", "mexia", 83),
                    new EdgeData("daytonaBeach", "jacksonville", 92),
                    new EdgeData("daytonaBeach", "orlando", 54),
                    new EdgeData("denver", "wichita", 523),
                    new EdgeData("denver", "grandJunction", 246),
                    new EdgeData("desMoines", "omaha", 135),
                    new EdgeData("desMoines", "minneapolis", 246),
                    new EdgeData("elPaso", "sanAntonio", 580),
                    new EdgeData("elPaso", "tucson", 320),
                    new EdgeData("eugene", "salem", 63),
                    new EdgeData("eugene", "medford", 165),
                    new EdgeData("europe", "philadelphia", 3939),
                    new EdgeData("ftWorth", "oklahomaCity", 209),
                    new EdgeData("fresno", "modesto", 109),
                    new EdgeData("grandJunction", "provo", 220),
                    new EdgeData("greenBay", "minneapolis", 304),
                    new EdgeData("greenBay", "milwaukee", 117),
                    new EdgeData("greensboro", "raleigh", 74),
                    new EdgeData("houston", "mexia", 165),
                    new EdgeData("indianapolis", "stLouis", 246),
                    new EdgeData("jacksonville", "savannah", 140),
                    new EdgeData("jacksonville", "lakeCity", 113),
                    new EdgeData("japan", "pointReyes", 5131),
                    new EdgeData("japan", "sanLuisObispo", 5451),
                    new EdgeData("kansasCity", "tulsa", 249),
                    new EdgeData("kansasCity", "stLouis", 256),
                    new EdgeData("kansasCity", "wichita", 190),
                    new EdgeData("keyWest", "tampa", 446),
                    new EdgeData("lakeCity", "tampa", 169),
                    new EdgeData("lakeCity", "tallahassee", 104),
                    new EdgeData("laredo", "sanAntonio", 154),
                    new EdgeData("laredo", "mexico", 741),
                    new EdgeData("lasVegas", "losAngeles", 275),
                    new EdgeData("lasVegas", "saltLakeCity", 486),
                    new EdgeData("lincoln", "wichita", 277),
                    new EdgeData("lincoln", "omaha", 58),
                    new EdgeData("littleRock", "memphis", 137),
                    new EdgeData("littleRock", "tulsa", 276),
                    new EdgeData("losAngeles", "sanDiego", 124),
                    new EdgeData("losAngeles", "sanLuisObispo", 182),
                    new EdgeData("medford", "redding", 150),
                    new EdgeData("memphis", "nashville", 210),
                    new EdgeData("miami", "westPalmBeach", 67),
                    new EdgeData("midland", "toledo", 82),
                    new EdgeData("minneapolis", "winnipeg", 463),
                    new EdgeData("modesto", "stockton", 29),
                    new EdgeData("montreal", "ottawa", 132),
                    new EdgeData("newHaven", "providence", 110),
                    new EdgeData("newHaven", "stamford", 92),
                    new EdgeData("newOrleans", "pensacola", 268),
                    new EdgeData("newYork", "philadelphia", 101),
                    new EdgeData("norfolk", "richmond", 92),
                    new EdgeData("norfolk", "raleigh", 174),
                    new EdgeData("oakland", "sanFrancisco", 8),
                    new EdgeData("oakland", "sanJose", 42),
                    new EdgeData("oklahomaCity", "tulsa", 105),
                    new EdgeData("orlando", "westPalmBeach", 168),
                    new EdgeData("orlando", "tampa", 84),
                    new EdgeData("ottawa", "toronto", 269),
                    new EdgeData("pensacola", "tallahassee", 120),
                    new EdgeData("philadelphia", "pittsburgh", 319),
                    new EdgeData("philadelphia", "newYork", 101),
                    new EdgeData("philadelphia", "uk1", 3548),
                    new EdgeData("philadelphia", "uk2", 3548),
                    new EdgeData("phoenix", "tucson", 117),
                    new EdgeData("phoenix", "yuma", 178),
                    new EdgeData("pointReyes", "redding", 215),
                    new EdgeData("pointReyes", "sacramento", 115),
                    new EdgeData("portland", "seattle", 174),
                    new EdgeData("portland", "salem", 47),
                    new EdgeData("reno", "saltLakeCity", 520),
                    new EdgeData("reno", "sacramento", 133),
                    new EdgeData("richmond", "washington", 105),
                    new EdgeData("sacramento", "sanFrancisco", 95),
                    new EdgeData("sacramento", "stockton", 51),
                    new EdgeData("salinas", "sanJose", 31),
                    new EdgeData("salinas", "sanLuisObispo", 137),
                    new EdgeData("sanDiego", "yuma", 172),
                    new EdgeData("saultSteMarie", "thunderBay", 442),
                    new EdgeData("saultSteMarie", "toronto", 436),
                    new EdgeData("seattle", "vancouver", 115),
                    new EdgeData("thunderBay", "winnipeg", 440)
            );

    private static class EdgeData {
        public String firstCityName;
        public String secondCityName;
        public int distance;

        public EdgeData(String firstCityName, String secondCityName, int distance) {
            this.firstCityName = firstCityName;
            this.secondCityName = secondCityName;
            this.distance = distance;
        }
    }

    private static class Vertex {
        public String cityName;
        public double latitude;
        public double longitude;

        public Vertex(String cityName, double lat, double longitude) {
            this.cityName = cityName;
            this.latitude = lat;
            this.longitude = longitude;
        }
    }

    private static class Edge {
        public Vertex firstCity;
        public Vertex secondCity;
        public int distance;

        public Edge(Vertex firstCity, Vertex secondCity, int distance) {
            this.firstCity = firstCity;
            this.secondCity = secondCity;
            this.distance = distance;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Edge edge = (Edge) o;
            return distance == edge.distance &&
                    Objects.equals(firstCity, edge.firstCity) &&
                    Objects.equals(secondCity, edge.secondCity);
        }

        @Override
        public int hashCode() {
            return Objects.hash(firstCity, secondCity, distance);
        }
    }

    private static class Graph {
        public Map<Vertex, Set<Edge>> adjList;
        public Map<String, Vertex> nameToVertex;
        public int minimumCostArc;

        public Graph() {
            this.adjList = new HashMap<>();
            this.nameToVertex = new HashMap<>();
            this.minimumCostArc = Integer.MAX_VALUE;
        }

        public Graph addVertex(Vertex v) {
            Set<Edge> edges = Optional.ofNullable(adjList.get(v)).orElseGet(HashSet::new);
            adjList.put(v, edges);
            nameToVertex.put(v.cityName, v);
            return this;
        }

        public Graph addEdge(Edge e) {
            Vertex v1 = e.firstCity;
            Vertex v2 = e.secondCity;
            Set<Edge> edges1 = Optional.ofNullable(adjList.get(v1)).orElseGet(HashSet::new);
            Set<Edge> edges2 = Optional.ofNullable(adjList.get(v2)).orElseGet(HashSet::new);
            edges1.add(e);
            Edge reverseEdge = new Edge(v2, v1, e.distance);
            edges2.add(reverseEdge);
            adjList.put(v1, edges1);
            adjList.put(v2, edges2);
            this.minimumCostArc = Math.min(minimumCostArc, e.distance);
            return this;
        }
    }

    public static Graph intializeGraph() {
        Graph graph = new Graph();
        for (Vertex v: VERTICES) {
            graph.addVertex(v);
        }
        for (EdgeData ed: EDGES) {
            Vertex firstV = graph.nameToVertex.get(ed.firstCityName);
            Vertex secondV = graph.nameToVertex.get(ed.secondCityName);
            graph.addEdge(new Edge(firstV, secondV, ed.distance));
        }
        return graph;
    }

    public static double sq(double x) {
        return x*x;
    }

    public static double sphericalHeuristic(Vertex v1, Vertex v2) {
        double lat1 = v1.latitude;
        double lat2 = v2.latitude;
        double long1 = v1.longitude;
        double long2 = v2.longitude;
        //sqrt(  (69.5 * (Lat1 - Lat2)  ) ^ 2 + (69.5 * cos((Lat1 + Lat2)/360 * pi) * (Long1 - Long2)) ^ 2).
        return Math.sqrt(sq(69.5 * (lat1 - lat2))
                + sq(69.5 * (Math.cos(((lat1+lat2)/360.0)*Math.PI)) * (long1 - long2)));
    }

    public static Map<Vertex, Long> getHopMap(Graph g, String source, String destination) {
        Vertex sourceV = g.nameToVertex.get(source);
        Vertex destV = g.nameToVertex.get(destination);
        Map<Vertex, Long> minimumDistance = g.nameToVertex.values().stream()
                .collect(Collectors.toMap(v -> v, v -> INF));
        Set<Vertex> discovered = new HashSet<>();
        minimumDistance.put(destV, 0L);
        Queue<Vertex> queue = new LinkedList<>();
        queue.add(destV);
        discovered.add(destV);
        while (!queue.isEmpty()) {
            Vertex current = queue.poll();
            List<Vertex> neighbors = g.adjList.get(current).stream().map(e -> e.secondCity)
                    .collect(Collectors.toList());
            for (Vertex neighbor: neighbors) {
                if (!discovered.contains(neighbor)) {
                    discovered.add(neighbor);
                    minimumDistance.put(neighbor, minimumDistance.get(current) + 1);
                    queue.add(neighbor);
                }
            }
        }
        return minimumDistance;
    }

    public static Stats DFS(Graph g, String source, String destination) {
        Vertex sourceV = g.nameToVertex.get(source);
        Vertex destV = g.nameToVertex.get(destination);
        Map<Vertex, Boolean> visited = new HashMap<>();
        for (Vertex v: g.adjList.keySet()) {
            visited.put(v, false);
        }
        if (sourceV.equals(destV)) {
            System.out.println("DFS");
            System.out.println("No. of cities expanded = " + 0);
            System.out.println("Max. length of queue during search = " + 1);
            System.out.println("Final path length = " + 0);
            System.out.println("Path = " + source);
            return new Stats(0, 1);
        }
        Reference<Long> maxQueueSize = new Reference<>(0L);
        List<Edge> path = DFSUtil(g, sourceV, destV, visited, 1, maxQueueSize);
        if (path == null) {
            System.err.println("Destination not reached");
            return null;
        }
        Collections.reverse(path);
        long noOfCitiesExpanded = visited.values().stream().filter(b -> b).count();
        long pathLength = path.stream().map(e -> e.distance).mapToLong(i -> i).sum();
        String pathString = Stream.concat(Stream.of(source), path.stream().map(e -> e.secondCity.cityName))
                .collect(Collectors.joining(","));
        System.out.println("DFS");
        System.out.println("No. of cities expanded = " + noOfCitiesExpanded);
        System.out.println("Max. length of queue during search = " + maxQueueSize.value);
        System.out.println("Final path length = " + pathLength);
        System.out.println("Path = " + pathString);
        return new Stats(noOfCitiesExpanded, maxQueueSize.value);
    }

    public static List<Edge> DFSUtil(Graph g, Vertex source, Vertex destination, Map<Vertex, Boolean> visited,
                                     long currentQueueSize, Reference<Long> maxQueueSize) {
        visited.put(source, true);
        Set<Edge> neighbors = g.adjList.get(source);
        Map<Vertex, Edge> neighborToEdge = neighbors.stream().collect(Collectors.toMap(e -> e.secondCity, e -> e));
        currentQueueSize += neighbors.size();
        maxQueueSize.update(Math.max(maxQueueSize.value, currentQueueSize));
        if (neighborToEdge.containsKey(destination)) {
            List<Edge> path = new ArrayList<>();
            path.add(neighborToEdge.get(destination));
            return path;
        }
        for (Edge e: neighbors) {
            if (!e.firstCity.equals(source)) {
                System.err.println("Error in edge direction or graph creation");
                return null;
            }
            Vertex otherV = e.secondCity;
            if (!visited.get(otherV)) {
                List<Edge> path = DFSUtil(g, otherV, destination, visited, currentQueueSize, maxQueueSize);
                if (path != null) {
                    path.add(e);
                    return path;
                }
            }
        }
        return null;
    }

    public static Stats AStar(Graph g, Function<Vertex, Double> hFun, String source, String destination) {
        Vertex sourceV = g.nameToVertex.get(source);
        Vertex destV = g.nameToVertex.get(destination);
        Set<Vertex> expanded = new HashSet<>();
        //Set<Vertex> inQueue = new HashSet<>();
        //Map<Vertex, Boolean> expanded = g.nameToVertex.values().stream().collect(Collectors.toMap(v -> v, v -> false));
        //Map<Vertex, Boolean> inQueue = g.nameToVertex.values().stream().collect(Collectors.toMap(v -> v, v -> false));
        Map<Vertex, Long> gFun = g.nameToVertex.values().stream()
                .collect(Collectors.toMap(v -> v, v -> INF));
        gFun.put(sourceV, 0L);

        Function<Vertex, Double> f = v -> gFun.get(v) + hFun.apply(v);
        PriorityQueue<Vertex> pq = new PriorityQueue<>(Comparator.comparing(f));
        pq.add(sourceV);
        int maxQueueSize = 1;
        int effectiveQueueSize = 1;
        //inQueue.add(sourceV);
        Map<Vertex, Vertex> previousVertexInShortestPath = new HashMap<>();
        previousVertexInShortestPath.put(sourceV, sourceV);
        while (!pq.isEmpty()) {
            Vertex head = pq.poll();
            if (expanded.contains(head)) {
                continue;
            } else {
                effectiveQueueSize--;
            }
            if (head.equals(destV)) {
                break;
            }
            expanded.add(head);
            Set<Edge> incidentEdges = g.adjList.get(head);
            long pathFromSourceToHead = gFun.get(head);
            for (Edge e: incidentEdges) {
                Vertex otherV = e.secondCity;
                long pathFromSourceToOtherV = gFun.get(otherV);
                long newPathLength = pathFromSourceToHead + e.distance;
                if (pathFromSourceToOtherV == INF) {
                    gFun.put(otherV, newPathLength);
                    pq.add(otherV);
                    effectiveQueueSize++;
                    maxQueueSize = Math.max(maxQueueSize, effectiveQueueSize);
                    previousVertexInShortestPath.put(otherV, head);
                } else if (pathFromSourceToOtherV > newPathLength) {
                    gFun.put(otherV, newPathLength);
                    previousVertexInShortestPath.put(otherV, head);
                    if (!expanded.contains(otherV)) {
                        pq.add(otherV);
                    }
                }
            }
        }
        List<Vertex> finalPath = new ArrayList<>();
        finalPath.add(destV);
        Vertex v = destV;
        while(true) {
            Vertex parent = previousVertexInShortestPath.get(v);
            if (parent == null || parent.equals(v)) {
                break;
            }
            finalPath.add(parent);
            v = parent;
        }
        Collections.reverse(finalPath);
        /*
        long temppl = 0L;
        for (int i = 0; i < finalPath.size() - 1; i++) {
            Vertex curr = finalPath.get(i);
            Vertex next = finalPath.get(i+1);
            Set<Edge> edges = g.adjList.get(curr);
            for (Edge e: edges) {
                if (e.secondCity.equals(next)) {
                    temppl += e.distance;
                    break;
                }
            }
        }
        System.out.println("temppl = " + temppl);
        */
        String pathString = finalPath.stream().map(x -> x.cityName).collect(Collectors.joining(","));
        System.out.println("A*");
        System.out.println("No. of cities expanded = " + expanded.size());
        System.out.println("Max. length of queue during search = " + maxQueueSize);
        System.out.println("Final path length = " + gFun.get(destV));
        System.out.println("Path = " + pathString);
        return new Stats(expanded.size(), maxQueueSize);
    }

    public static Stats RBFS(Graph g, Function<Vertex, Double> hFun, String source, String destination) {
        final long INF = 100000000000000L;
        Vertex sourceV = g.nameToVertex.get(source);
        Vertex destV = g.nameToVertex.get(destination);
        Reference<Long> noOfCitiesExpanded = new Reference<>(0L);
        Reference<Long> maxQueueSize = new Reference<>(0L);
        Node result = RBFSUtil(g, new Node(sourceV, null, 0, hFun.apply(sourceV)), null,
                destV, INF, hFun, noOfCitiesExpanded, 0, maxQueueSize);
        if (result == null) {
            System.err.println("Destination not found by RBFS");
            return null;
        }

        List<Vertex> pathList = new ArrayList<>();
        pathList.add(result.vertex);
        Node current = result;
        while(true) {
            Node parent = current.parent;
            if (parent == null) {
                break;
            }
            pathList.add(parent.vertex);
            current = parent;
        }
        Collections.reverse(pathList);
        String pathString = pathList.stream().map(v -> v.cityName).collect(Collectors.joining(","));

        System.out.println("RBFS");
        System.out.println("No. of cities expanded = " + noOfCitiesExpanded.value);
        System.out.println("Max. length of queue during search = " + maxQueueSize.value);
        System.out.println("Final path length = " + result.gValue);
        System.out.println("Path = " + pathString);
        return new Stats(noOfCitiesExpanded.value, maxQueueSize.value);
    }

    public static Node RBFSUtil(Graph g, Node current, Node parent,
                                Vertex destination,
                                double fLimit, Function<Vertex,
                                Double> hFun, Reference<Long> expandCount,
                                long successorQueueSize, Reference<Long> maxQueueSize) {
        if (current.vertex.equals(destination)) {
            return current;
        }
        expandCount.update(expandCount.value + 1);
        List<Node> children = current.generateChildren(g, parent, hFun);
        successorQueueSize += children.size();
        maxQueueSize.update(Math.max(maxQueueSize.value, successorQueueSize));
        if (children.isEmpty()) {
            current.fValue = INF;
            return null;
        }
        for (Node child: children) {
            child.fValue = Math.max(child.gValue + child.hValue, current.fValue);
        }
        while (true) {
            children.sort(Node::compareTo);
            Node bestNode = children.get(0);
//            System.out.println("Best node = " + bestNode.vertex.cityName);
//            System.out.println("Best f-value = " + bestNode.fValue);
            if (bestNode.fValue > fLimit) {
                current.fValue = bestNode.fValue;
                return null;
            }
            Node result;
            if (children.size() > 1) {
                result = RBFSUtil(g, bestNode, current, destination, Math.min(fLimit, children.get(1).fValue),
                        hFun, expandCount, successorQueueSize, maxQueueSize);
            } else {
                result = RBFSUtil(g, bestNode, current, destination, fLimit, hFun, expandCount, successorQueueSize,
                        maxQueueSize);
            }
            if (result != null) {
                return result;
            }
        }
    }

    private static class Node implements Comparable<Node> {
        public Vertex vertex;
        public Node parent;
        public double fValue;
        public long gValue;
        public double hValue;

        public Node(Vertex v, Node parent, long gValue, double hValue) {
            this.vertex = v;
            this.parent = parent;
            this.gValue = gValue;
            this.hValue = hValue;
            this.fValue = gValue + hValue;
        }

        public List<Node> generateChildren(Graph g, Node parent, Function<Vertex, Double> hFun) {
            return g.adjList.get(vertex).stream()
                                        .filter(e -> (parent == null || !e.secondCity.equals(parent.vertex)))
                                        .map(e -> new Node(e.secondCity, this, gValue + e.distance,
                                                hFun.apply(e.secondCity)))
                                        .collect(Collectors.toList());
        }

        @Override
        public int compareTo(Node node) {
            return Double.compare(fValue, node.fValue);
        }
    }

    public static void main(String[] args) {
        Graph g = intializeGraph();
        String algo = args[0];
        String heuristicPref = args[1];
        String source = args[2];
        String dest = args[3];

        if (algo.equals("DFS")) {
            DFS(g, source, dest);
        } else if (algo.equals("A*")) {
            if (heuristicPref.equals("0")) {
                AStar(g, vertex -> sphericalHeuristic(vertex, g.nameToVertex.get(dest)), source, dest);
            } else {
                Map<Vertex, Long> hops = getHopMap(g, source, dest);
                Function<Vertex, Double> hopHeuristic =
                        vertex -> (((double) g.minimumCostArc) * hops.get(vertex));
                AStar(g, hopHeuristic, source, dest);
            }
        } else if (algo.equals("RBFS")) {
            if (heuristicPref.equals("0")) {
                RBFS(g, vertex -> sphericalHeuristic(vertex, g.nameToVertex.get(dest)), source, dest);
            } else {
                Map<Vertex, Long> hops = getHopMap(g, source, dest);
                Function<Vertex, Double> hopHeuristic =
                        vertex -> (((double) g.minimumCostArc) * hops.get(vertex));
                RBFS(g, hopHeuristic, source, dest);
            }
        } else {
            System.out.println("Unrecognized algo = " + algo);
        }
    }

    public static void analyze(Graph g) {
        List<Vertex> vertices = new ArrayList<>(g.nameToVertex.values());
        int noOfVertices = vertices.size();

        for (int firstIndex = 0; firstIndex < noOfVertices; firstIndex++) {
            for (int secondIndex = 0; secondIndex < noOfVertices; secondIndex++) {
                if (firstIndex != secondIndex) {
                    Vertex dest = vertices.get(secondIndex);
                    Vertex source = vertices.get(firstIndex);
                    Stats stats1 = AStar(g, v -> sphericalHeuristic(v, dest), vertices.get(firstIndex).cityName,
                            dest.cityName);
                    Map<Vertex, Long> hops = getHopMap(g, source.cityName, dest.cityName);
                    Function<Vertex, Double> hopHeuristic =
                            vertex -> (((double) g.minimumCostArc) * hops.get(vertex));
                    Stats stats2 = AStar(g, hopHeuristic, source.cityName, dest.cityName);
                    if (stats1.expandedCities != stats2.expandedCities) {
                        System.out.println(String.format("Expanded cities are not equal for (%s,%s): (%d, %d)",
                                source.cityName, dest.cityName, stats1.expandedCities, stats2.expandedCities));
                    }
                    if (stats1.maxQueueSize != stats2.maxQueueSize) {
                        System.out.println(String.format("Max queue sizes cities are not equal for (%s,%s): (%d, %d)",
                                source.cityName, dest.cityName, stats1.maxQueueSize, stats2.maxQueueSize));
                    }
                }
            }
        }

        /*
        long expandedCitySum = 0L;
        long expandedCityMax = Long.MIN_VALUE;
        long expandedCityMin = Long.MAX_VALUE;
        long queueSizeSum = 0L;
        long queueSizeMax = Long.MIN_VALUE;
        long queueSizeMin = Long.MAX_VALUE;
        Pair<Integer,Integer> worstPairForExpandedCities = null;
        Pair<Integer,Integer> worstPairForQueueSize = null;
        long count = 0;
        for (int firstIndex = 0; firstIndex < noOfVertices; firstIndex++) {
            for (int secondIndex = 0; secondIndex < noOfVertices; secondIndex++) {
                if (firstIndex != secondIndex) {
                    Stats stats = DFS(g, vertices.get(firstIndex).cityName, vertices.get(secondIndex).cityName);
                    if (stats.expandedCities > expandedCityMax) {
                        expandedCityMax = stats.expandedCities;
                        worstPairForExpandedCities = Pair.of(firstIndex, secondIndex);
                    }
                    expandedCityMin = Math.min(expandedCityMin, stats.expandedCities);
                    expandedCitySum += stats.expandedCities;
                    if (stats.maxQueueSize > queueSizeMax) {
                        queueSizeMax = stats.maxQueueSize;
                        worstPairForQueueSize = Pair.of(firstIndex, secondIndex);
                    }
                    queueSizeMin = Math.min(queueSizeMin, stats.maxQueueSize);
                    queueSizeSum += stats.maxQueueSize;
                    count++;
                }
            }
        }
        System.out.println("For DFS");
        System.out.println("Number of cities expanded during a single path search (Max) = " + expandedCityMax);
        System.out.println("Number of cities expanded during a single path search (Min) = " + expandedCityMin);
        System.out.println("Number of cities expanded during a single path search (Avg) = " + expandedCitySum / count);
        System.out.println("Maximum size of the queue during a single path search (Max) = " + queueSizeMax);
        System.out.println("Maximum size of the queue during a single path search (Min) = " + queueSizeMin);
        System.out.println("Maximum size of the queue during a single path search (Avg) = " + queueSizeSum / count);
        System.out.println(String.format("Worst pair for expanded cities = (%s,%s)",
                vertices.get(worstPairForExpandedCities.getLeft()).cityName,
                vertices.get(worstPairForExpandedCities.getRight()).cityName));
        System.out.println(String.format("Worst pair for queue size = (%s,%s)",
                vertices.get(worstPairForQueueSize.getLeft()).cityName,
                vertices.get(worstPairForQueueSize.getRight()).cityName));


        expandedCitySum = 0L; //TODO: This may overflow
        expandedCityMax = Long.MIN_VALUE;
        expandedCityMin = Long.MAX_VALUE;
        queueSizeSum = 0L;
        queueSizeMax = Long.MIN_VALUE;
        queueSizeMin = Long.MAX_VALUE;
        count = 0;
        vertices = new ArrayList<>(g.nameToVertex.values());
        noOfVertices = vertices.size();
        for (int firstIndex = 0; firstIndex < noOfVertices; firstIndex++) {
            for (int secondIndex = 0; secondIndex < noOfVertices; secondIndex++) {
                if (firstIndex != secondIndex) {
                    Vertex dest = vertices.get(secondIndex);
                    Stats stats = AStar(g, v -> sphericalHeuristic(v, dest), vertices.get(firstIndex).cityName,
                            dest.cityName);
                    if (stats.expandedCities > expandedCityMax) {
                        expandedCityMax = stats.expandedCities;
                        worstPairForExpandedCities = Pair.of(firstIndex, secondIndex);
                    }
                    expandedCityMin = Math.min(expandedCityMin, stats.expandedCities);
                    expandedCitySum += stats.expandedCities;
                    if (stats.maxQueueSize > queueSizeMax) {
                        queueSizeMax = stats.maxQueueSize;
                        worstPairForQueueSize = Pair.of(firstIndex, secondIndex);
                    }
                    queueSizeMin = Math.min(queueSizeMin, stats.maxQueueSize);
                    queueSizeSum += stats.maxQueueSize;
                    count++;
                }
            }
        }
        System.out.println("For A*");
        System.out.println("Number of cities expanded during a single path search (Max) = " + expandedCityMax);
        System.out.println("Number of cities expanded during a single path search (Min) = " + expandedCityMin);
        System.out.println("Number of cities expanded during a single path search (Avg) = " + expandedCitySum / count);
        System.out.println("Maximum size of the queue during a single path search (Max) = " + queueSizeMax);
        System.out.println("Maximum size of the queue during a single path search (Min) = " + queueSizeMin);
        System.out.println("Maximum size of the queue during a single path search (Avg) = " + queueSizeSum / count);
        System.out.println(String.format("Worst pair for expanded cities = (%s,%s)",
                vertices.get(worstPairForExpandedCities.getLeft()).cityName,
                vertices.get(worstPairForExpandedCities.getRight()).cityName));
        System.out.println(String.format("Worst pair for queue size = (%s,%s)",
                vertices.get(worstPairForQueueSize.getLeft()).cityName,
                vertices.get(worstPairForQueueSize.getRight()).cityName));


        expandedCitySum = 0L; //TODO: This may overflow
        expandedCityMax = Long.MIN_VALUE;
        expandedCityMin = Long.MAX_VALUE;
        queueSizeSum = 0L;
        queueSizeMax = Long.MIN_VALUE;
        queueSizeMin = Long.MAX_VALUE;
        count = 0;
        vertices = new ArrayList<>(g.nameToVertex.values());
        noOfVertices = vertices.size();
        for (int firstIndex = 0; firstIndex < noOfVertices; firstIndex++) {
            for (int secondIndex = 0; secondIndex < noOfVertices; secondIndex++) {
                if (firstIndex != secondIndex) {
                    Vertex dest = vertices.get(secondIndex);
                    Stats stats = RBFS(g, v -> sphericalHeuristic(v, dest), vertices.get(firstIndex).cityName,
                            dest.cityName);
                    if (stats.expandedCities > expandedCityMax) {
                        expandedCityMax = stats.expandedCities;
                        worstPairForExpandedCities = Pair.of(firstIndex, secondIndex);
                    }
                    expandedCityMin = Math.min(expandedCityMin, stats.expandedCities);
                    expandedCitySum += stats.expandedCities;
                    if (stats.maxQueueSize > queueSizeMax) {
                        queueSizeMax = stats.maxQueueSize;
                        worstPairForQueueSize = Pair.of(firstIndex, secondIndex);
                    }
                    queueSizeMin = Math.min(queueSizeMin, stats.maxQueueSize);
                    queueSizeSum += stats.maxQueueSize;
                    count++;
                }
            }
        }
        System.out.println("For RBFS");
        System.out.println("Number of cities expanded during a single path search (Max) = " + expandedCityMax);
        System.out.println("Number of cities expanded during a single path search (Min) = " + expandedCityMin);
        System.out.println("Number of cities expanded during a single path search (Avg) = " + expandedCitySum / count);
        System.out.println("Maximum size of the queue during a single path search (Max) = " + queueSizeMax);
        System.out.println("Maximum size of the queue during a single path search (Min) = " + queueSizeMin);
        System.out.println("Maximum size of the queue during a single path search (Avg) = " + queueSizeSum / count);
        System.out.println(String.format("Worst pair for expanded cities = (%s,%s)",
                vertices.get(worstPairForExpandedCities.getLeft()).cityName,
                vertices.get(worstPairForExpandedCities.getRight()).cityName));
        System.out.println(String.format("Worst pair for queue size = (%s,%s)",
                vertices.get(worstPairForQueueSize.getLeft()).cityName,
                vertices.get(worstPairForQueueSize.getRight()).cityName));
        */
    }

    private static class Stats {
        public long expandedCities;
        public long maxQueueSize;

        public Stats(long expandedCities, long maxQueueSize) {
            this.expandedCities = expandedCities;
            this.maxQueueSize = maxQueueSize;
        }
    }

    public static class Reference<T> {
        public T value;

        public Reference(T value) {
            this.value = value;
        }

        public void update(T val) {
            value = val;
        }
    }

    public static class Pair<L, R> {
        private L left;
        private R right;

        private Pair(L left, R right) {
            this.left = left;
            this.right = right;
        }

        public L getLeft() {
            return left;
        }

        public R getRight() {
            return right;
        }

        public static <A, B> Pair<A, B> of(A a, B b) {
            return new Pair<>(a, b);
        }

        @Override
        public String toString() {
            return "Pair{" +
                    "left=" + left +
                    ", right=" + right +
                    '}';
        }
    }
}
