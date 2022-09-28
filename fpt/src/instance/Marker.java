package instance;

import java.util.*;

public class Marker {
    private static int nextId = 0;

    Integer id;
    private Integer chromosome;
    Integer sign;
    Integer startPosition;
    Integer prev_ir;
    Integer next_ir;
    Marker nextFixed;
    public final String name;
    public final Integer value;
    private Set<Marker> matches;
    Set<Integer> matchIds;

    Marker(int id, int sign, int startPosition, int name, Genome g, int prev_ir, int next_ir) {
        this.id = id;
        this.chromosome = 1;
        this.name = Integer.toString(name);
        this.value = name;
        this.startPosition = startPosition;
        this.sign = sign;
        this.matches = new HashSet<Marker>();
        this.matchIds = new HashSet<Integer>();
        this.prev_ir = prev_ir;
        this.next_ir = next_ir;
        g.idGeneMap.put(id, this);
    }

    Marker(Marker oldMarker, Genome g) {
        this.id = oldMarker.id;
        this.chromosome = oldMarker.chromosome;
        this.name = oldMarker.name;
        this.value = oldMarker.value;
        this.startPosition = oldMarker.startPosition;
        this.sign = oldMarker.sign;
        this.prev_ir = oldMarker.prev_ir;
        this.next_ir = oldMarker.next_ir;
        this.matches = new HashSet<Marker>();
        this.matchIds = new HashSet<Integer>();
        for (Marker m : oldMarker.getMatches()) {
            this.matchIds.add(m.id);
        }
        g.idGeneMap.put(this.id, this);
    }

    @Override
    public boolean equals(Object o) {
        if(o instanceof Marker){
            Marker toCompare = (Marker) o;
            return this.id.equals(toCompare.id);
          }
          return false;
    }

    @Override
    public int hashCode() {
        return this.id;
    }

    public boolean isFixed() {
        if (getMatches().size() == 1 && this.getMatch().getMatches().size() == 1) {
            return true;
        } else {
            return false;
        }
    }

    public void print() {
        System.out.println("Tag = " + this.tag());
        // System.out.println("Accession Number = "+ this.name);
        // System.out.println("Chromosome = "+ this.getChromosome());
        // System.out.println("Sign = "+ this.sign);
        // System.out.println("Start = "+ this.startPosition);
        if (this.isFixed()) {
            System.out.println("Fixed with " + this.getMatch().tag());
        } else {
            System.out.println("Not fixed with possibilities:");
            for (Iterator<Integer> it = matchIds.iterator(); it.hasNext();) {
                System.out.println("Possible: " + it.next());
            }
        }

    }

    public Integer getId() {
        return this.id.intValue();
    }

    public void addMatch(Marker g) {
        this.matchIds.add(g.id);
        this.getMatches().add(g);
    }

    public Marker getMatch() {
        if (this.getMatches().isEmpty())
            return null;
        else
            return this.getMatches().iterator().next();
    }

    public Integer getValue() {
        return this.value * this.sign.intValue();
    }

    public Integer getNextIR() {
        return this.next_ir;
    }

    public Integer getPrevIR() {
        return this.prev_ir;
    }

    public boolean isForwardMatch(Marker g, Boolean first) {
        if (this.getMatches().contains(g) && this.sign.intValue() == g.sign.intValue()
                && (first || this.prev_ir.intValue() == g.prev_ir.intValue()))
            return true;
        else
            return false;
    }

    public boolean isForwardMatchLeft(Marker g) {
        if (this.getMatches().contains(g) && this.sign.intValue() == g.sign.intValue()
                && this.next_ir.intValue() == g.next_ir.intValue())
            return true;
        else
            return false;
    }

    public boolean isBackwardMatch(Marker g, Boolean first) {
        if (this.getMatches().contains(g) && this.sign.intValue() != g.sign.intValue()
                && (first || this.prev_ir.intValue() == g.next_ir.intValue()))
            return true;
        else
            return false;
    }

    public boolean isBackwardMatchRight(Marker g) {
        if (this.getMatches().contains(g) && this.sign.intValue() != g.sign.intValue()
                && this.next_ir.intValue() == g.prev_ir.intValue())
            return true;
        else
            return false;
    }

    public boolean isMatch(Marker m2) {
        return this.getMatches().contains(m2);
    }

    public int matchNum() {
        return this.getMatches().size();
    }

    public int getChromosome() {
        return chromosome.intValue();
    }

    public void setChromosome(Integer chromosome) {
        this.chromosome = chromosome;
    }

    public String tag() {
        return this.name + "." + this.id;
    }

    public Set<Marker> getMatches() {
        return matches;
    }

    public void setMatches(Set<Marker> matches) {
        this.matches = matches;
    }

    public boolean isRare() {
        Marker match = this.getMatch();
        if (match == null)
            return false;
        else
            return (match.matchNum() <= this.matchNum());
    }

    public boolean isRare(int del) {
        Marker match = this.getMatch();
        if (match == null)
            return false;
        else
            return (match.matchNum() <= this.matchNum() - del);
    }

}
