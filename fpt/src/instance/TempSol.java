package instance;

import instance.FixedMatch.Added;

import java.util.*;
import main.MCSP;

public class TempSol {
    Set<FixedMatch> matches;
    int k;

    public class ReducedStr {
        public String str;
        public int cost;
    }

    public TempSol() {
        k = 0;
        this.matches = new HashSet<FixedMatch>();
    }

    /*
     * public void addFix(Marker m1, Marker m2){
     * this.matches.add(new FixedMatch(m1, m2));
     * }
     * 
     * public void addFix(FixedMatch fm){
     * this.matches.add(fm);
     * }
     */
    public void addFix(Marker m1, Marker m2, boolean incrementBlockCount) {
        if (incrementBlockCount)
            k++;
        this.matches.add(new FixedMatch(m1, m2));
    }

    public void addFix(Marker m1, Marker m2, boolean incrementBlockCount, Added where) {
        if (incrementBlockCount)
            k++;
        this.matches.add(new FixedMatch(m1, m2, where));
    }

    public void addFix(FixedMatch fm, boolean incrementBlockCount) {
        if (incrementBlockCount)
            k++;
        this.matches.add(fm);
    }

    public TempSol copyTempSol() {
        TempSol newTemp = new TempSol();
        for (FixedMatch fm : matches) {
            newTemp.addFix(fm, false);
        }
        newTemp.k = k;
        return newTemp;
    }

    public int fixNumber() {
        return this.matches.size();
    }

    public void print() {
        for (FixedMatch fm : matches) {
            System.out.println(fm.m1.tag() + "--" + fm.m2.tag() + " " + fm.where);
        }
    }

    public void checkIsComplete(Instance inst) throws Exception {

        /*
         * To be tested.
         * Remark: k can be bigger than the actual number of blocks, that's fine.
         * 
         */
        final int unmatched = -2;
        int n1 = inst.getGenome1().size();
        int n2 = inst.getGenome2().size();
        int[] matchedTo1 = new int[n1];
        boolean[] matched2 = new boolean[n2];

        // init
        for (int i = 0; i < n1; i++) {
            matchedTo1[i] = unmatched;
        }
        for (int i = 0; i < n2; i++) {
            matched2[i] = false;
        }
        // check all matches
        for (FixedMatch fm : matches) {
            Marker mark1 = inst.genome1.getMarkerById(fm.m1.getId());
            Marker mark2 = inst.genome2.getMarkerById(fm.m2.getId());
            if (mark1 == null || mark2 == null) {
                mark1 = inst.genome1.getMarkerById(fm.m2.getId());
                mark2 = inst.genome2.getMarkerById(fm.m1.getId());
            }
            if (mark1 == null || mark2 == null)
                throw new Exception("matched markers from the same genome or problem with Ids, fm : " + fm.toString());
            int pos1 = inst.genome1.getMarkerIndex(mark1);
            int pos2 = inst.genome2.getMarkerIndex(mark2);

            if (matchedTo1[pos1] != unmatched)
                throw new Exception("double match in 1 at pos" + pos1);
            if (matched2[pos2])
                throw new Exception("double match in 2 at pos " + pos2);
            matchedTo1[pos1] = pos2;
            matched2[pos2] = true;
        }
        int prec = unmatched;
        int computedK = 0;
        for (int i = 0; i < n1; i++) { // signed case: change this loop
            int curr = matchedTo1[i];
            // compute CSP size: count non deleted markers which are not parallel to prec
            if (curr != unmatched && curr != prec + 1)
                computedK++;
            // rare markers in genome 1
            if (curr == unmatched && inst.getGenome1().getMarker(i).isRare())
                throw new Exception("unmatched rare marker in genome 1, pos=" + i);
            prec = curr;
        }
        // rare markers in genome 2
        for (int i = 0; i < n2; i++) {
            if (!matched2[i] && inst.getGenome2().getMarker(i).isRare())
                throw new Exception("unmatched rare marker in genome 2, pos=" + i);
        }

        if (computedK > k)
            throw new Exception("wrong number of blocks");
        /**/

    }

    public void incrementBlockCount() {
        k++;
    }

    public int getK() {
        return k;
    }

    public ReducedStr getReducedGenomes(Genome g1, Genome g2) throws Exception {
        ReducedStr ans = new ReducedStr();
        ans.cost = 0;
        Integer v = 1;
        final int unmatched = -10;

        int n1 = g1.size();
        int n2 = g2.size();
        Block[] blocks1 = new Block[n1];
        Block[] blocks2 = new Block[n2];
        int[] irs1 = new int[n1];
        int[] irs2 = new int[n2];
        int[] matchedTo1 = new int[n1];
        boolean[] matched2 = new boolean[n2];

        for (int i = 0; i < n1; i++) {
            matchedTo1[i] = unmatched;
            irs1[i] = unmatched;
            blocks1[i] = null;
        }
        for (int i = 0; i < n2; i++) {
            matched2[i] = false;
            irs2[i] = unmatched;
            blocks2[i] = null;
        }

        for (FixedMatch fm : matches) {
            // System.out.println(fm);
            Marker mark1 = g1.getMarkerById(fm.m1.getId());
            Marker mark2 = g2.getMarkerById(fm.m2.getId());
            if (mark1 == null || mark2 == null) {
                mark1 = g1.getMarkerById(fm.m2.getId());
                mark2 = g2.getMarkerById(fm.m1.getId());
            }
            if (mark1 == null || mark2 == null)
                throw new Exception("matched markers from the same genome or problem with Ids, fm : " + fm.toString());

            int pos1 = g1.getMarkerIndex(mark1);
            int pos2 = g2.getMarkerIndex(mark2);

            if (matchedTo1[pos1] != unmatched)
                throw new Exception("double match in 1 at pos" + pos1);
            if (matched2[pos2])
                throw new Exception("double match in 2 at pos " + pos2);
            matchedTo1[pos1] = pos2;
            matched2[pos2] = true;
        }
        int prec = unmatched;
        int computedK = 0;
        int pos1 = -1;
        int pos2 = -1;
        boolean rev = false;
        for (int i = 0; i < n1; i++) {
            int curr = matchedTo1[i];
            // compute CSP size: count non deleted markers which are not parallel to prec
            if (curr == unmatched) { // current value on deleted block
                if (prec == unmatched) { /// predecessor already on deleted block
                    blocks1[pos1].addGene(g1.getMarker(i).getValue());
                    blocks1[pos1].addIR(g1.getMarker(i).getPrevIR());
                } else { // predecessor on common block
                    pos1 = i;
                    blocks1[pos1] = new Block(true);
                    blocks1[pos1].addGene(g1.getMarker(i).getValue());
                    irs1[pos1] = g1.getMarker(i).getPrevIR();
                    if (rev) {
                        irs2[pos2] = g2.getMarker(prec).getPrevIR();
                    }
                }
            } else if (!rev && curr == prec + 1 && g1.getMarker(i).isForwardMatch(g2.getMarker(curr), false) || rev && curr == prec - 1 && g1.getMarker(i).isBackwardMatch(g2.getMarker(curr), false)) { // current in the same common block as previous
                blocks1[pos1].addGene(g1.getMarker(i).getValue());
                blocks1[pos1].addIR(g1.getMarker(i).getPrevIR());
                blocks2[pos2].addGene(g2.getMarker(curr).getValue());
                if (rev) {
                    blocks2[pos2].addIR(g2.getMarker(curr).getNextIR());
                } else {
                    blocks2[pos2].addIR(g2.getMarker(curr).getPrevIR());
                }
            } else { // current in another common block
                if (rev && prec != unmatched) {
                    irs2[pos2] = g2.getMarker(prec).getPrevIR();
                }
                pos1 = i;
                pos2 = curr;
                int a1 = g1.getMarker(pos1).getValue();
                int a2 = g2.getMarker(pos2).getValue();
                blocks1[pos1] = new Block(false);
                blocks1[pos1].addGene(a1);
                irs1[pos1] = g1.getMarker(pos1).getPrevIR();
                if (a1 == a2) {
                    blocks2[pos2] = new Block(false, false);
                    blocks2[pos2].addGene(a2);
                    irs2[pos2] = g2.getMarker(pos2).getPrevIR();
                    rev = false;
                } else {
                    blocks2[pos2] = new Block(false, true);
                    blocks2[pos2].addGene(a2);
                    rev = true;
                }
                computedK++;
            }
            // rare markers in genome 1
            if (curr == unmatched && g1.getMarker(i).isRare()) {
                MCSP.printV("unmatched rare marker in genome 1, pos=" + i);
                return null;
            }
            prec = curr;
        }
        // rare markers in genome 2
        Boolean prev_not_matched = false;
        for (int i = 0; i < n2; i++) {
            if (!matched2[i]) {
                if (prev_not_matched) {
                    blocks2[pos2].addGene(g2.getMarker(i).getValue());
                    blocks2[pos2].addIR(g2.getMarker(i).getPrevIR());
                } else {
                    prev_not_matched = true;
                    pos2 = i;
                    blocks2[pos2] = new Block(true);
                    blocks2[pos2].addGene(g2.getMarker(i).getValue());
                    irs2[pos2] = g2.getMarker(i).getPrevIR();
                }
                if (g2.getMarker(i).isRare()) {
                    MCSP.printV("unmatched rare marker in genome 2, pos=" + i);
                    return null;
                }
            } else {
                prev_not_matched = false;
            }
        }

        String sg1 = "", sg2 = "", si1 = "", si2 = "";

        if (MCSP.reduce) {
            HashMap<Block.InvBlock, Integer> blockMap = new HashMap<Block.InvBlock, Integer>();
            Boolean first = true;
            for (int i = 0; i < n1; i++) {
                if (blocks1[i] != null) {
                    Block.InvBlock b = new Block.InvBlock(blocks1[i]);
                    Integer v_;
                    if (b.isDeleted()) {
                        v_ = v;
                        v++;
                    } else {
                        v_ = blockMap.putIfAbsent(b, v);
                        if (v_ == null) {
                            v_ = v;
                            v++;
                        }
                    }
                    if (first) {
                        first = false;
                    } else {
                        MCSP.print(Integer.toString(irs1[i]));
                        si1 += " " + Integer.toString(irs1[i]);
                    }
                    b.print();
                    if (b.isNegative()) {
                        sg1 += " " + Integer.toString(-1 * v_);
                    } else {
                        sg1 += " " + Integer.toString(v_);
                    }
                    ans.cost++;
                }
            }

            MCSP.print("--------------------------");

            first = true;
            for (int i = 0; i < n2; i++) {
                if (blocks2[i] != null) {
                    Block.InvBlock b = new Block.InvBlock(blocks2[i]);
                    Integer v_;
                    if (b.isDeleted()) {
                        v_ = v;
                        v++;
                    } else {
                        v_ = blockMap.putIfAbsent(b, v);
                        if (v_ == null) {
                            v_ = v;
                            v++;
                        }
                    }
                    if (first) {
                        first = false;
                    } else {
                        MCSP.print(Integer.toString(irs2[i]));
                        si2 += " " + Integer.toString(irs2[i]);
                    }
                    b.print();
                    if (b.isNegative()) {
                        sg2 += " " + Integer.toString(-1 * v_);
                    } else {
                        sg2 += " " + Integer.toString(v_);
                    }
                    ans.cost++;
                }
            }
        } else {
            HashMap<Block.InvBlock, ArrayList<Integer>> blockMap = new HashMap<Block.InvBlock, ArrayList<Integer>>();
            Boolean first = true;
            for (int i = 0; i < n1; i++) {
                if (blocks1[i] != null) {
                    Block.InvBlock b = new Block.InvBlock(blocks1[i]);
                    Integer v_;
                    if (b.isDeleted()) {
                        v_ = v;
                        v++;
                    } else {
                        blockMap.putIfAbsent(b, new ArrayList<Integer>());
                        ArrayList<Integer> l = blockMap.get(b);
                        l.add(v);
                        v_ = v;
                        v += b.size();
                    }
                    if (first) {
                        first = false;
                    } else {
                        MCSP.print(Integer.toString(irs1[i]));
                        si1 += " " + Integer.toString(irs1[i]);
                    }
                    b.print();

                    ArrayList<Integer> genes = new ArrayList<Integer>();
                    for (int j = 0; j < b.size(); j++) {
                        genes.add(v_);
                        v_++;
                    }
                    if (b.isNegative()) {
                        Collections.reverse(genes);
                    }
                    for (int j = 0; j < b.size() - 1; j++) {
                        sg1 += " " + Integer.toString(genes.get(j));
                        si1 += " " + Integer.toString(b.getIR(j));
                    }
                    sg1 += " " + Integer.toString(genes.get(b.size() - 1));
                    ans.cost++;
                }
            }

            MCSP.print("--------------------------");

            first = true;
            for (int i = 0; i < n2; i++) {
                if (blocks2[i] != null) {
                    Block.InvBlock b = new Block.InvBlock(blocks2[i]);
                    Integer v_;
                    if (b.isDeleted()) {
                        v_ = v;
                        v++;
                    } else {
                        ArrayList<Integer> l = blockMap.getOrDefault(b, new ArrayList<Integer>());
                        if (l.isEmpty()) {
                            v_ = v;
                            v += b.size();
                        } else {
                            v_ = l.get(0);
                            l.remove(0);
                        }
                    }
                    if (first) {
                        first = false;
                    } else {
                        MCSP.print(Integer.toString(irs2[i]));
                        si2 += " " + Integer.toString(irs2[i]);
                    }
                    b.print();

                    ArrayList<Integer> genes = new ArrayList<Integer>();
                    for (int j = 0; j < b.size(); j++) {
                        genes.add(v_);
                        v_++;
                    }
                    if (b.isNegative()) {
                        Collections.reverse(genes);
                    }
                    for (int j = 0; j < b.size() - 1; j++) {
                        sg2 += " " + Integer.toString(genes.get(j));
                        si2 += " " + Integer.toString(b.getIR(j));
                    }
                    sg2 += " " + Integer.toString(genes.get(b.size() - 1));
                    ans.cost++;
                }
            }
        }

        if (MCSP.commonCost) {
            ans.cost = computedK;
        } else {
            ans.cost = ans.cost - 1 - computedK;
        }
        ans.str = sg1 + "\n" + si1 + "\n" + sg2 + "\n" + si2;
        return ans;
    }

}
