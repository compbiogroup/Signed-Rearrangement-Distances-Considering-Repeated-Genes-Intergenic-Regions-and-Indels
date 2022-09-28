package algorithm;

import java.util.*;

import main.MCSP;
import instance.*;
import instance.FixedMatch.Added;

public class ReductionRules {
    Instance inst;
    Genome genome1, genome2;
    int matchReduced = 0;
    public static boolean ApplyRule1 = true;
    public static boolean ApplyRule2 = true;
    public static boolean ApplyRule3 = true;
    public static boolean ApplyRule4 = true;

    public ReductionRules(Instance inst) {
        this.inst = inst;
        this.genome1 = inst.getGenome1();
        this.genome2 = inst.getGenome2();
    }

    // returns the number of removed matches because of P3 Rule
    public int applyRules() throws Exception {
        updateGenome();
        matchReduced = 0;
        boolean reducible = true;
        while (reducible) {
            reducible = false;
            /*
             * boolean success= true;
             * if (ApplyRule1) {
             * while (success){
             * success = BorderRule();
             * MCSP.printV("Border Success");
             * if (success) reducible=true;
             * }
             * success= true;
             * }
             * if (ApplyRule2) {
             * while (success){
             * success = ParallelRule();
             * MCSP.printV("Parallel Success");
             * if (success) reducible=true;
             * }
             * success= true;
             * }
             * if (ApplyRule3) {
             * while (success){
             * success = reduceP3Matches(genome1);
             * MCSP.printV("P3 Success");
             * if (success) {
             * reducible=true;
             * j++;
             * }
             * }
             * success= true;
             * }
             */
            if (ApplyRule1) {
                if (applyRule(1))
                    reducible = true;
            }
            if (ApplyRule2) {
                if (applyRule(2))
                    reducible = true;
            }
            if (ApplyRule3) {
                if (applyRule(3))
                    reducible = true;
            }
            if (ApplyRule4) {
                if (applyRule(4))
                    reducible = true;
            }
        }

        // inst.printFull();
        return matchReduced;
    }

    private boolean applyRule(int i) throws Exception {
        boolean reducible = false;
        boolean success = true;
        while (success) {
            success = false;
            switch (i) {
                case 1:
                    success = BorderRule();
                    break;
                case 2:
                    success = ParallelRule();
                    break;
                case 3: {
                    success = reduceStarMatches(genome1);
                    if (!success)
                        success = reduceStarMatches(genome2);
                    if (success)
                        matchReduced++;
                    break;
                }
                case 4: {
                    success = reduceK22Matches();
                    if (success)
                        matchReduced += 2;
                    break;
                }
                default:
                    break;
            }
            if (success)
                reducible = true;
        }
        return reducible;
    }

    public boolean BorderRule() {
        boolean hit = false;
        Marker currentMatch = null;
        Marker previousMatch = null;
        Marker nextMatch = null;
        Marker current = null;
        Marker previous = null;
        Marker next = null;
        ListIterator<Marker> it = genome1.iterator();
        if (it.hasNext())
            previous = it.next();
        if (it.hasNext())
            current = it.next();
        if (it.hasNext())
            next = it.next();
        if (next == null)
            return false;
        while (next != null && !hit) {
            MCSP.printV("Current = " + current.name);
            if (current.isFixed() && previous.isFixed() && next.isFixed()) {
                currentMatch = current.getMatch();
                previousMatch = genome2.getPredecessor(currentMatch);
                nextMatch = genome2.getSuccessor(currentMatch);
                if (previousMatch != null && nextMatch != null) {
                    if (currentMatch.isFixed() &&
                            previousMatch.isFixed() &&
                            nextMatch.isFixed()// &&
                    // !next.isMatch(previousMatch)&&
                    // !previous.isMatch(nextMatch)&&
                    // !previous.isMatch(genome2.getPredecessor(next.getMatch()))&&
                    // !previousMatch.isMatch(genome1.getPredecessor(nextMatch.getMatch()))
                    )
                        hit = true;
                }
            }
            if (it.hasNext() && !hit) {
                previous = current;
                current = next;
                next = it.next();
            } else if (!hit)
                return false;
        }
        if (hit) {
            MCSP.printV("Trying to remove marker: " + current.tag());
            if (!inst.checkParallel(previous, current) &&
                    !inst.checkParallel(current, next)) {
                foundNewBlock();
                removeAndFix(current, currentMatch, Added.DR1);
                MCSP.printV("Is a block");
                if (inst.checkParallel(previous, next))
                    foundNewBlock();
                if (inst.checkParallel(previousMatch.getMatch(), nextMatch.getMatch()))
                    foundNewBlock();
            } else {
                MCSP.printV("Removing because parallel to previous or next");
                removeAndFix(current, currentMatch, Added.DR2);
            }
        }
        return hit;
    }

    public boolean reduceStarMatches(Genome g) {
        boolean hit = false;
        for (ListIterator<Marker> it = g.iterator(); (it.hasNext() && !hit);) {
            Marker m = it.next();
            if (m.matchNum() == 1) {
                Iterator<Marker> itMatches = m.getMatches().iterator();
                Marker m1 = itMatches.next();
                // Marker m2 = itMatches.next();
                if (m1.matchNum() > 1 && matchGivesTwoBreakpoints(m, m1, g)) {
                    m.getMatches().remove(m1);
                    m1.getMatches().remove(m);
                    hit = true;
                }
                // else if (matchGivesTwoBreakpoints(m,m2,g)){
                // m.getMatches().remove(m2);
                // m2.getMatches().remove(m);
                // hit=true;
                // }
                // }
            }
        }
        if (hit)
            MCSP.printV("Successfully Applied P3 Rule");
        return hit;
    }

    private boolean matchGivesTwoBreakpoints(Marker middle, Marker end, Genome middleGenome) {
        Marker pMiddle, sMiddle, pEnd, sEnd;
        Genome endGenome = inst.getOther(middleGenome);
        pMiddle = middleGenome.getPredecessor(middle);
        sMiddle = middleGenome.getSuccessor(middle);
        pEnd = endGenome.getPredecessor(end);
        sEnd = endGenome.getSuccessor(end);
        if (pMiddle == null || pEnd == null || sMiddle == null || sEnd == null)
            return false;
        if (middle.isForwardMatch(end, true)) {
            if (!pMiddle.isForwardMatchLeft(pEnd) && !sMiddle.isForwardMatch(sEnd, false)) {
                return true;
            }
        } else {
            if (!pMiddle.isBackwardMatchRight(sEnd) && !sMiddle.isBackwardMatch(pEnd, false)) {
                return true;
            }
        }
        return false;
    }

    public boolean reduceK22Matches() {
        boolean hit = false;
        for (ListIterator<Marker> it = genome1.iterator(); (it.hasNext() && !hit);) {
            Marker m1 = it.next();
            if (m1.matchNum() == 2 && m1.getMatch().matchNum() == 2) {
                Iterator<Marker> itMatches = m1.getMatches().iterator();
                Marker n1 = itMatches.next();
                Marker n2 = itMatches.next();
                Marker m2 = null;
                for (Marker u : n1.getMatches()) {
                    if (u != m1)
                        m2 = u;
                }
                if (matchGivesTwoBreakpoints(m1, n1, genome1) && matchGivesTwoBreakpoints(m2, n2, genome1)) {
                    m1.getMatches().remove(n1);
                    n1.getMatches().remove(m1);
                    m2.getMatches().remove(n2);
                    n2.getMatches().remove(m2);
                    return true;
                }
            }
        }
        return hit;
    }

    // returns true in case the rule has been successfully applied
    public boolean ParallelRule() throws Exception {
        boolean hit = false;
        int indexG1;
        Marker g = null;
        // findNextFixed
        for (ListIterator<Marker> it = genome1.iterator(); (it.hasNext() && !hit);) {
            indexG1 = it.nextIndex();
            g = it.next();
            Marker g2 = null;
            if (g.isFixed() && it.hasNext()) {
                // find next fixed gene
                boolean reachedNextFixed = false;
                ListIterator<Marker> it2 = genome1.iterator(indexG1 + 1);
                while (!reachedNextFixed && it2.hasNext()) {
                    g2 = it2.next();
                    if (g2.isFixed()) {
                        reachedNextFixed = true;
                    }
                }
                if (reachedNextFixed)
                    if (inst.checkParallel(g, g2))
                        hit = true;

            }
        }

        if (hit) {
            Marker succ = genome1.getSuccessor(g);
            Marker succM;
            if (g.isForwardMatch(g.getMatch(), true))
                succM = genome2.getSuccessor(g.getMatch());
            else if (g.isBackwardMatch(g.getMatch(), true))
                succM = genome2.getPredecessor(g.getMatch());
            else
                throw new Exception("No match when in hit of parallel rule");
            removeAndFix(succ, succM, Added.DR3);
        }
        return hit;
    }

    private void foundNewBlock() {
        this.inst.getTempSol().incrementBlockCount();
        // inst.foundBlock();
    }

    // TODO Reduction rules for monotone blocks

    private void updateGenome() {
        this.genome1 = inst.getGenome1();
        this.genome2 = inst.getGenome2();
    }

    private void removeAndFix(Marker m1, Marker m2, Added where) {
        genome1.remove(m1);
        genome2.remove(m2);
        this.inst.getTempSol().addFix(m1, m2, false, where);
    }

    public static void switchRule(boolean apply, int ruleNum) {
        switch (ruleNum) {
            case 1:
                ApplyRule1 = apply;
                break;
            case 2:
                ApplyRule2 = apply;
                break;
            case 3:
                ApplyRule3 = apply;
                break;
            case 4:
                ApplyRule4 = apply;
                break;
            default:
                break;
        }
    }
}
