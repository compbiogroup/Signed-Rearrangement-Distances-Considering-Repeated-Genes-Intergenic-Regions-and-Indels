package algorithm;

import utilities.Timer;
import instance.*;
import instance.FixedMatch.Added;

import java.util.*;

import main.MCSP;

public class SampleGraph {
    Instance inst;
    boolean isNoInstance;
    Set<Marker> V1, V2;
    LinkedHashSet<Marker> V;
    Set<Edge> blackEdges;
    List<Marker> blackVertices;
    HashMap<Marker, Set<SampleGraph.Edge>> incidentEdges;
    HashMap<Marker, Integer> deletions;

    enum Color {
        black, green, red, blue, orange
    };

    class Edge {
        Color color_up;
        Color color_down;
        Marker v1, v2;
        boolean forward;

        Edge(Color c1, Color c2, Marker v1, Marker v2, boolean f) {
            this.color_up = c1;
            this.color_down = c2;
            this.v1 = v1;
            this.v2 = v2;
            this.forward = f;
        }

        public String toString() {
            return (" " + v1.tag() + "-" + v2.tag() + " " + color_up + " " + color_down);
        }
    }

    class DelMarker {
        Marker m;

        DelMarker() {
            m = null;
        }
    }

    class Extremes {
        private Set<Marker> seen;
        private HashMap<Marker, Marker> pairs;
        private HashMap<Marker, DelMarker> deletion_spots;

        Extremes() {
            seen = new HashSet<Marker>();
            pairs = new HashMap<Marker, Marker>();
            deletion_spots = new HashMap<Marker, DelMarker>();
        }

        public Set<Marker> markers() {
            return new HashSet<Marker>(seen);
        }

        public void setDeleteSpot(Marker m) {
            deletion_spots.get(m).m = m;
        }

        public Marker getDeleteSpot(Marker m) {
            return deletion_spots.get(m).m;
        }

        public boolean inDeletionZone(Marker m) {
            return deletion_spots.containsKey(m);
        }

        public Marker getPair(Marker m) {
            return pairs.get(m);
        }

        public boolean isSeen(Marker m) {
            return seen.contains(m) || seen.contains(getPair(m));
        }

        public void remove(Marker m) {
            Marker pm = getPair(m);
            seen.remove(pm);
            seen.remove(m);
        }

        public void restore(Marker m) {
            if (!isSeen(m)) {
                seen.add(m);
            }
        }

        public boolean isEmpty() {
            return seen.isEmpty();
        }

        public void print() {
            for (Marker m : seen) {
                System.out.print(m.tag() + " ");
            }
            System.out.println();
        }

        public void make_spots(Genome g1, Genome g2) throws Exception {
            for (Marker m1_ : seen) {
                Marker ms[] = {m1_, getPair(m1_)};
                for (Marker m1 : ms) {
                    if (!deletion_spots.containsKey(m1)) {
                        Genome g;
                        if (g1.contains(m1)) {
                            g = g1;
                        } else if (g2.contains(m1)) {
                            g = g2;
                        } else {
                            throw new Exception("Every marker has to be in one of the genomes in make_spots.");
                        }
                        ListIterator<Marker> it = g.iterator(g.getMarkerIndex(m1));
                        while (it.hasPrevious()) {
                            Marker m2 = it.previous();
                            if (isVertex(m2) && !this.isSeen(m2) && degree(m2) != 0) {
                                it.next();
                                break;
                            }
                        }
                        DelMarker dm = new DelMarker();
                        while (it.hasNext()) {
                            Marker m2 = it.next();
                            if (isVertex(m2) && !this.isSeen(m2) && degree(m2) != 0) {
                                break;
                            }
                            deletion_spots.put(m2, dm);
                        }
                    }
                }
            }
        }
    }

    public SampleGraph() {
        this.V1 = new HashSet<Marker>();
        this.V2 = new HashSet<Marker>();
        this.V = new LinkedHashSet<Marker>();
        this.blackEdges = new HashSet<SampleGraph.Edge>();
        this.blackVertices = new LinkedList<Marker>();
        this.incidentEdges = new HashMap<Marker, Set<SampleGraph.Edge>>();
        this.deletions = new HashMap<Marker, Integer>();
        this.isNoInstance = false;
    }

    public SampleGraph(SampleGraph sg) {
        inst = InstanceFactory.copyInstance(sg.inst);
        isNoInstance = sg.isNoInstance;
        V1 = sg.V1;
        V2 = sg.V2;
        V = sg.V;
        blackEdges = new HashSet<Edge>(sg.blackEdges);
        blackVertices = new LinkedList<Marker>(sg.blackVertices);
        incidentEdges = new HashMap<Marker, Set<SampleGraph.Edge>>();
        for (Marker m : sg.incidentEdges.keySet()) {
            incidentEdges.put(m, new HashSet<SampleGraph.Edge>(sg.incidentEdges.get(m)));
        }
        deletions = new HashMap<Marker, Integer>(sg.deletions);
    }

    public boolean hasBlackParallel() {
        boolean parallel = false;
        ListIterator<Marker> it = blackVertices.listIterator();
        while (it.hasNext() && !parallel) {
            Marker m = it.next();
            ListIterator<Marker> it2 = blackVertices.listIterator(blackVertices.indexOf(m));
            while (it2.hasNext() && !parallel) {
                Marker m2 = it2.next();
                if (inst.checkParallel(m, m2))
                    parallel = true;
            }
        }
        return parallel;
    }

    SampleGraph(Instance inst, TempSol.ReducedStr best_sol) throws Exception {
        this();
        MCSP.print("Building Sample Graph");
        this.inst = inst;
        Genome genome1 = inst.getGenome1();
        Genome genome2 = inst.getGenome2();
        // initialize vertex sets
        addGenomeMarkers(genome1, V1);
        addGenomeMarkers(genome2, V2);
        // initialize edges
        for (Iterator<Marker> it = V.iterator(); it.hasNext();) {
            incidentEdges.put(it.next(), new HashSet<SampleGraph.Edge>());
        }
        // add one color at a time
        MCSP.print("Adding black edges");
        addBlackEdges(best_sol);
        addGreenEdges(genome1, genome2);
        addRedEdges(genome1, genome2);
    }

    private void addGreenEdges(Genome genome1, Genome genome2) {
        MCSP.printV("Trying green for black edge");
        for (Iterator<Edge> it = blackEdges.iterator(); it.hasNext();) {
            Edge be = it.next();
            if (be.forward) { // green for forward matches
                MCSP.printV("Forward black edge " + be.toString());
                Iterator<Marker> mIt = genome1.startFromMarker(be.v1);
                Iterator<Marker> mIt2 = genome2.startFromMarker(be.v2);
                Marker point = mIt.next();
                Marker point2 = mIt2.next();
                while (mIt.hasNext() && mIt2.hasNext() && point.isForwardMatchLeft(point2)) {
                    point = mIt.next();
                    point2 = mIt2.next();
                    if (point.isForwardMatch(point2, true) && !(blackVertices.contains(point))
                            && !(blackVertices.contains(point2))) {
                        Edge ge = new Edge(Color.green, Color.green, point, point2, true);
                        addEdge(point, point2, ge);
                        MCSP.printV("Added green edge" + ge.toString());
                    }
                }
            } else { // green for backward matches
                MCSP.printV("Backward black edge " + be.toString());
                ListIterator<Marker> mIt = genome1.startFromMarker(be.v1);
                ListIterator<Marker> mIt2 = genome2.startFromMarker(be.v2);
                Marker point = mIt.next();
                Marker point2 = mIt2.next();
                point2 = mIt2.previous();
                while (mIt.hasNext() && mIt2.hasPrevious() && point.isBackwardMatchRight(point2)) {
                    point = mIt.next();
                    point2 = mIt2.previous();
                    if (point.isBackwardMatch(point2, true) && !(blackVertices.contains(point))
                            && !(blackVertices.contains(point2))) {
                        Edge ge = new Edge(Color.green, Color.red, point, point2, false);
                        addEdge(point, point2, ge);
                        MCSP.printV("Added green edge" + ge.toString());
                    }
                }
            }
        }
    }

    private void addRedEdges(Genome genome1, Genome genome2) {
        MCSP.print("Adding red edges for black edge");
        for (Iterator<Edge> it = blackEdges.iterator(); it.hasNext();) {
            Edge be = it.next();
            if (be.forward) { // red only for forward matches
                ListIterator<Marker> mIt = genome1.startFromMarker(be.v1);
                ListIterator<Marker> mIt2 = genome2.startFromMarker(be.v2);
                Marker point = be.v1;
                Marker point2 = be.v2;
                while (mIt.hasPrevious() && mIt2.hasPrevious() && point.isForwardMatch(point2, false)) {
                    point = mIt.previous();
                    point2 = mIt2.previous();
                    if (point.isForwardMatchLeft(point2) && !(blackVertices.contains(point))
                            && !(blackVertices.contains(point2))) {
                        Edge re = new Edge(Color.red, Color.red, point, point2, true);
                        addEdge(point, point2, re);
                        MCSP.print("Added red edge " + re.toString());
                    }
                }
            } else {
                ListIterator<Marker> mIt = genome1.startFromMarker(be.v1);
                ListIterator<Marker> mIt2 = genome2.startFromMarker(be.v2);
                Marker point = be.v1;
                Marker point2 = mIt2.next();
                while (mIt.hasPrevious() && mIt2.hasNext() && point.isBackwardMatch(point2, false)) {
                    point = mIt.previous();
                    point2 = mIt2.next();
                    if (point.isBackwardMatchRight(point2) && !(blackVertices.contains(point))
                            && !(blackVertices.contains(point2))) {
                        Edge re = new Edge(Color.red, Color.green, point, point2, false);
                        addEdge(point, point2, re);
                        MCSP.print("Added red edge " + re.toString());
                    }
                }
            }
        }
    }

    private void addBlackEdges(TempSol.ReducedStr best_sol) throws Exception {
        // one endpoint must be in V1
        for (Iterator<Marker> it = V1.iterator(); it.hasNext() && !isNoInstance;) {
            Marker m = it.next();
            if (m.isFixed()) {
                Marker m2 = m.getMatch();
                if (m2.isFixed()) {
                    Edge be = new Edge(Color.black, Color.black, m, m2, m.isForwardMatch(m2, true));
                    addEdge(m, m2, be);
                    MCSP.printV("Added black edge" + be.toString());
                    blackEdges.add(be);
                    if (best_sol != null && blackEdges.size() - (MCSP.commonCost ? 0 : 1) >= best_sol.cost)
                        isNoInstance = true;
                    blackVertices.add(m);
                    blackVertices.add(m2);
                } else {
                    throw new Exception("Both markers must be fixed");
                }
            }
        }
        MCSP.print("Number of black edges: " + blackEdges.size());
        // inst.print();
    }

    private void addGenomeMarkers(Genome g, Set<Marker> s) {
        for (ListIterator<Marker> it = g.iterator(); it.hasNext();) {
            Marker m = it.next();
            s.add(m);
            V.add(m);
        }
    }

    private void addEdge(Marker v1, Marker v2, Edge e) {
        this.incidentEdges.get(v1).add(e);
        this.incidentEdges.get(v2).add(e);
    }

    public boolean isVertex(Marker m) {
        return V.contains(m);
    }

    public int degree(Marker m) {
        if (this.incidentEdges.containsKey(m))
            return this.incidentEdges.get(m).size();
        else
            return 0;
    }

    public void print() {
        System.out.println("-----------------------");
        for (Map.Entry<Marker, Set<SampleGraph.Edge>> e : incidentEdges.entrySet()) {
            System.out.print(e.getKey().tag() + " " + e.getKey().getNextIR() + ":");
            for (SampleGraph.Edge edge : e.getValue()) {
                System.out.print(edge.toString());
            }
            System.out.println();
        }
    }

    public Marker findIsolated() {
        Marker minDeg = null;
        for (Iterator<Marker> it = V.iterator(); it.hasNext();) {
            Marker m = it.next();
            if (degree(m) == 0 && m.isRare()) {
                if (minDeg == null) {
                    minDeg = m;
                } else if (degree(m) == 0 && minDeg.matchNum() > m.matchNum() && m.isRare()) {
                    minDeg = m;
                }
            }
        }
        if (minDeg != null)
            MCSP.print("Minimum Match Number of isolated: " + minDeg.matchNum());
        else
            MCSP.print("No isolated");
        return minDeg;
    }

    private Set<Marker> getNeighbors(Marker v) {
        if (!isVertex(v) || !incidentEdges.containsKey(v))
            return null;
        else {
            Set<Marker> neighbors = new HashSet<Marker>();
            for (Edge e : incidentEdges.get(v)) {
                neighbors.add(e.v1);
                neighbors.add(e.v2);
            }
            neighbors.remove(v);
            return neighbors;
        }
    }

    public Set<Marker> findOddPath(boolean rare) {
        if (rare) {
            MCSP.print("Looking for rare odd paths");
        } else {
            MCSP.print("Looking for odd paths");
        }
        Set<Marker> candidatePath = null;
        for (Marker v : V) {
            // find path starting on v, v must be rare
            if ((!rare || v.isRare()) && degree(v) == 1 && !blackVertices.contains(v)) {
                MCSP.print("Starting path at " + v.tag());
                Set<Marker> currentPath = new HashSet<Marker>();
                currentPath.add(v);
                Set<Marker> neighbors = getNeighbors(v);
                for (Marker m : neighbors) {
                    if (!currentPath.contains(m)) {// find growing end
                        v = m;
                        MCSP.print("Growing path with " + v.tag());
                    }
                }
                currentPath.addAll(neighbors);
                boolean mayGrow = true;
                boolean isCycle = false;
                while (mayGrow) {
                    if (degree(v) == 1) {
                        MCSP.print("Reached End Of Path");
                        mayGrow = false;
                    } else {
                        neighbors = getNeighbors(v);
                        boolean foundNew = false;
                        for (Marker m : neighbors) {
                            if (!currentPath.contains(m)) {// find growing end
                                v = m;
                                foundNew = true;
                                MCSP.print("Growing path with " + v.name);
                            }
                        }
                        if (!foundNew) {
                            isCycle = true;
                            MCSP.print("Is a cycle");
                            mayGrow = false;
                        }
                        currentPath.addAll(neighbors);
                    }
                }
                MCSP.print("Parity of path: " + (currentPath.size() % 2));
                if (!isCycle && !((currentPath.size() % 2) == 0)) {
                    if (candidatePath == null ||
                            candidatePath.size() * candidatePath.iterator().next().matchNum() > currentPath.size()
                                    * v.matchNum())
                        candidatePath = currentPath;
                    MCSP.print("Found odd path");
                }
            }
        }

        MCSP.print("Have not found odd path");
        return candidatePath;

    }

    public void addFixes(Timer timer, Boolean use_timer) throws Exception {
        for (Iterator<Edge> it = blackEdges.iterator(); it.hasNext();) {
            Edge edge = it.next();
            fixEdge(edge);
            it.remove();
        }
        // checkForDegreeOne
        boolean hasDegreeOne = true;
        while (hasDegreeOne) {
            hasDegreeOne = false;
            for (Marker v : V) {
                if (use_timer && timer.done()) {
                    throw new InterruptedException();
                }
                if (degree(v) == 1 && !blackVertices.contains(v)) {
                    hasDegreeOne = true;
                    // degree one and rare=>even path either red or green
                    if (v.isRare()) {
                        Edge edge = getEdge(v);
                        fixAndIsolateEdge(edge);
                    }
                    // degree one and abundant=>odd path take green edge
                    else {
                        Edge edge = getEdge(v);
                        if (inst.getGenome1().contains(v) && edge.color_up == Color.green ||
                            inst.getGenome2().contains(v) && edge.color_down == Color.green) {
                            fixAndIsolateEdge(edge);
                        }
                    }
                    MCSP.printV("Found Degree One");
                }
            }
        }
        MCSP.printV("Removed all paths");
        for (Marker v : V) {
            if (use_timer && timer.done()) {
                throw new InterruptedException();
            }
            if (degree(v) == 2) {
                Iterator<Edge> it = incidentEdges.get(v).iterator();
                Edge edge1 = it.next();
                Edge edge2 = it.next();
                if (edge1.color_up == Color.green) {
                    fixAndIsolateEdge(edge1);
                    MCSP.printV("recursing");
                    // recurse to find paths or other degree two
                    addFixes(timer, use_timer);
                } else if (edge2.color_up == Color.green) {
                    fixAndIsolateEdge(edge2);
                    // recurse to find paths or other degree two
                    MCSP.printV("recursing");
                    addFixes(timer, use_timer);
                } else {
                    throw new Exception("A vertex with two non-green edges should not exist");
                }
            }
        }
    }

    private Edge getEdge(Marker m) {
        if (this.incidentEdges.containsKey(m))
            return this.incidentEdges.get(m).iterator().next();
        else
            return null;
    }

    private void fixEdge(Edge e) {
        inst.getTempSol().addFix(e.v1, e.v2, e.color_up == Color.black, Added.End);
    }

    private void fixAndIsolateEdge(Edge e) {
        inst.getTempSol().addFix(e.v1, e.v2, e.color_up == Color.black, Added.End);
        makeIsolated(e.v1);
        makeIsolated(e.v2);
    }

    class EdgeData {
        Set<Edge> edges;
        HashMap<Marker, Set<SampleGraph.Edge>> edge_sets;
        Set<Marker> vertices;
        HashMap<Marker, Integer> deletions;

        EdgeData() {
            edges = new HashSet<Edge>();
            edge_sets = new HashMap<Marker, Set<SampleGraph.Edge>>();
            vertices = new HashSet<Marker>();
            deletions = new HashMap<Marker, Integer>();
        }
    }

    public boolean breakOddPath(Marker m, EdgeData ed) {
        Genome genome1 = inst.getGenome1();
        for (Edge e : incidentEdges.get(m)) {
            ListIterator<Marker> mIt = genome1.startFromMarker(e.v1);
            if (e.color_up == Color.green) {
                mIt.next();
                while (mIt.hasNext()) {
                    Marker v = mIt.next();
                    HashSet<Edge> greens = new HashSet<Edge>();
                    for (Edge eg : incidentEdges.get(v)) {
                        if (eg.color_up == Color.green) {
                            Marker[] vauxs = new Marker[2];
                            vauxs[0] = eg.v1;
                            vauxs[1] = eg.v2;
                            for (Marker vaux : vauxs) {
                                if (degree(vaux) == 1) {
                                    if (vaux.isRare(deletions.getOrDefault(vaux, 0))) {
                                        return false;
                                    } else {
                                        deletions.put(vaux, deletions.getOrDefault(vaux, 0) + 1);
                                        ed.deletions.put(vaux, ed.deletions.getOrDefault(vaux, 0) + 1);
                                    }
                                }
                            }
                            greens.add(eg);
                        }
                    }
                    if (greens.isEmpty()) {
                        break;
                    }
                    for (Edge eg : greens) {
                        incidentEdges.get(eg.v1).remove(eg);
                        incidentEdges.get(eg.v2).remove(eg);
                        ed.edges.add(eg);
                    }
                }
            } else if (e.color_up == Color.red) {
                while (mIt.hasPrevious()) {
                    Marker v = mIt.previous();
                    HashSet<Edge> reds = new HashSet<Edge>();
                    for (Edge er : incidentEdges.get(v)) {
                        if (er.color_up == Color.red) {
                            Marker[] vauxs = new Marker[2];
                            vauxs[0] = er.v1;
                            vauxs[1] = er.v2;
                            for (Marker vaux : vauxs) {
                                if (degree(vaux) == 1) {
                                    if (vaux.isRare(deletions.getOrDefault(vaux, 0))) {
                                        return false;
                                    } else {
                                        deletions.put(vaux, deletions.getOrDefault(vaux, 0) + 1);
                                        ed.deletions.put(vaux, ed.deletions.getOrDefault(vaux, 0) + 1);
                                    }
                                }
                            }
                            reds.add(er);
                        }
                    }
                    if (reds.isEmpty()) {
                        break;
                    }
                    for (Edge eg : reds) {
                        incidentEdges.get(eg.v1).remove(eg);
                        incidentEdges.get(eg.v2).remove(eg);
                        ed.edges.add(eg);
                    }
                }
            }
        }
        this.removeVertex(m,ed);
        return true;
    }

    public void restoreOddPath(EdgeData ed) {
        for (Edge e : ed.edges) {
            incidentEdges.get(e.v1).add(e);
            incidentEdges.get(e.v2).add(e);
        }
        for (Map.Entry<Marker, Integer> pair : ed.deletions.entrySet()) {
            deletions.put(pair.getKey(), deletions.get(pair.getKey()) - pair.getValue());
        }
        ed.deletions.clear();
        restoreVertices(ed);
    }


    public void removeVertex(Marker m, EdgeData ed) {
        Set<Marker> Nm = getNeighbors(m);
        Set<SampleGraph.Edge> edge_set = incidentEdges.get(m);
        incidentEdges.remove(m);
        if (Nm != null) {
            for (Marker n : Nm) {
                for (Iterator<Edge> it = incidentEdges.get(n).iterator(); it.hasNext();) {
                    Edge e = it.next();
                    if (e.v1 == m || e.v2 == m) {
                        it.remove();
                    }
                }
            }
        }
        ed.edge_sets.put(m, edge_set);
        ed.vertices.add(m);
    }

    public void restoreVertices(EdgeData ed) {
        for (Marker m : ed.vertices) {
            incidentEdges.put(m, ed.edge_sets.get(m));
            for (Edge e : ed.edge_sets.get(m)) {
                if (e.v1 == m) {
                    this.incidentEdges.get(e.v2).add(e);
                }
                if (e.v2 == m) {
                    this.incidentEdges.get(e.v1).add(e);
                }
            }
        }
        ed.vertices.clear();
        ed.edge_sets.clear();
    }

    private void makeIsolated(Marker m) {
        Set<Marker> Nm = getNeighbors(m);
        incidentEdges.remove(m);
        for (Marker n : Nm) {
            for (Iterator<Edge> it = incidentEdges.get(n).iterator(); it.hasNext();) {
                Edge e = it.next();
                if (e.v1 == m || e.v2 == m)
                    it.remove();
            }
        }
    }

}
