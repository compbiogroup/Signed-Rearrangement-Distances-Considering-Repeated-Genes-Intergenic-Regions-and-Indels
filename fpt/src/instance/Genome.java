package instance;

import java.util.*;

import main.MCSP;

public class Genome {
    ArrayList<Marker> markers;
    Map<Integer, Marker> idGeneMap;
    TempSol currentSol;
    int id;

    public Genome(int id) {
        this.id = id;
        markers = new ArrayList<Marker>();
        idGeneMap = new HashMap<Integer, Marker>();
    }

    public void addMarker(Marker g) {
        this.markers.add(g);
        this.idGeneMap.put(g.id, g);
    }

    public boolean contains(Marker m) {
        if (this.markers.contains(m))
            return true;
        else
            return false;
    }

    public Marker getMarker(int i) {
        if (i > this.markers.size() - 1 || i < 0) {
            return null;
        }
        return this.markers.get(i);
    }

    public Marker getMarkerById(Integer i) {
        return this.idGeneMap.get(i);
    }

    public ListIterator<Marker> iterator() {
        return this.markers.listIterator();
    }

    public ListIterator<Marker> iterator(int index) {
        return this.markers.listIterator(index);
    }

    public ListIterator<Marker> startFromMarker(Marker m) {
        return this.markers.listIterator(this.getMarkerIndex(m));
    }

    public List<Marker> geneList() {
        return this.markers;
    }

    public void testGenome(Genome other) {
        for (ListIterator<Marker> it = markers.listIterator(); it.hasNext();) {
            Marker g = it.next();
            MCSP.printV("Checking Marker at " + g.startPosition);
            if (it.hasNext()) {
                Marker g2 = markers.get(it.nextIndex());
                // System.out.println("G2 "+ g2.startPosition);
                if (g2.startPosition < g.startPosition && (g2.getChromosome() <= g.getChromosome())) {
                    System.out.println("We have a problem with");
                    System.out.println("Gene 1 starts at " + g.startPosition);
                    System.out.println("Gene 2 starts at " + g2.startPosition);
                }
            }
        }
        checkMatches(other);
    }

    public void removeNonMatched() {
        int i = 0;
        for (ListIterator<Marker> it = markers.listIterator(); it.hasNext();) {
            Marker g = it.next();
            if (g.getMatches().size() == 0) {
                it.remove();
                i++;
            }
        }
        MCSP.print("Removed " + i + " non-Matched Markers");
    }

    public void removeNegative() {
        MCSP.print("Removing Negative from one Genome");
        Set<Marker> negMarker = new HashSet<Marker>();

        for (Iterator<Marker> it = this.iterator(); it.hasNext();) {
            Marker m = it.next();
            if (m.sign.intValue() <= 0) {
                negMarker.add(m);
            }
        }
        for (Marker m : negMarker) {
            this.remove(m);
        }
    }

    public void checkMatches(Genome other) {
        for (Iterator<Marker> it = this.markers.iterator(); it.hasNext();) {
            Marker g = it.next();
            Marker match = g.getMatch();
            if (match != null && other.getMarkerIndex(match) < 0) {
                System.out.println("Alarmmmm");
            }
        }
    }

    public int getMarkerIndex(Marker g) {
        return this.markers.indexOf(g);
    }

    public Marker getSuccessor(Marker g) {
        return this.getMarker(this.getMarkerIndex(g) + 1);
    }

    public Marker getPredecessor(Marker g) {
        return this.getMarker(this.getMarkerIndex(g) - 1);
    }

    public void remove(Marker g) {
        ListIterator<Marker> git = this.iterator(this.getMarkerIndex(g));
        Marker prev = null;
        Marker next = null;
        int ir = 0;
        if (git.hasPrevious()) {
            prev = git.previous();
            ir += prev.next_ir;
            git.next();
        }
        git.next();
        if (git.hasNext()) {
            next = git.next();
            ir += next.prev_ir;
            next.prev_ir = ir;
        }
        if (prev != null) {
            prev.next_ir = ir;
        }
        this.markers.remove(g);
        for (Iterator<Marker> it = g.getMatches().iterator(); it.hasNext();) {
            Marker match = it.next();
            match.getMatches().remove(g);
        }
    }

    public void info() {
        System.out.println("Genome " + id + "; Size: " + markers.size());
        System.out.println("DMax: " + dMax() + "; DAvg: " + dAvg());
    }

    public void print() {
        this.info();
        String markerString = "";
        for (Iterator<Marker> it = markers.iterator(); it.hasNext();) {
            Marker m = it.next();
            // m.print();
            String sign;
            if (m.sign > 0) {
                sign = "";
            } else
                sign = "-";
            markerString += " (" + sign + m.tag() + ") " + m.getNextIR();
        }
        System.out.println(markerString);
    }

    public int dMax() {
        int dMax = -1;
        for (ListIterator<Marker> it = markers.listIterator(); it.hasNext();) {
            Marker m = it.next();
            dMax = Math.max(dMax, m.matchNum());
        }
        return dMax;
    }

    public int size() {
        return this.markers.size();
    }

    public float dAvg() {
        int dSum = 0;
        for (ListIterator<Marker> it = markers.listIterator(); it.hasNext();) {
            Marker m = it.next();
            dSum += m.matchNum();
        }
        return (float) dSum / (float) markers.size();
    }

    public void makePointersTransitive(Genome other) {
        boolean applicable = true;
        Marker m1, m2, m3; // m1 and m3 have degree one, m1 and m2 are matches of m
        while (applicable) {
            applicable = false;
            for (ListIterator<Marker> it = markers.listIterator(); it.hasNext();) {
                Marker m = it.next();
                if (m.matchNum() > 1) {
                    for (Iterator<Marker> matchIt = m.getMatches().iterator(); matchIt.hasNext();) {
                        m1 = matchIt.next();
                        for (Iterator<Marker> matchIt2 = m.getMatches().iterator(); matchIt2.hasNext();) {
                            m2 = matchIt2.next();
                            if (m2.matchNum() > 1) {
                                for (Iterator<Marker> matchIt3 = m2.getMatches().iterator(); matchIt3.hasNext();) {
                                    m3 = matchIt3.next();
                                    if (!m1.isMatch(m3)) {
                                        m1.addMatch(m3);
                                        m3.addMatch(m1);
                                        applicable = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}
