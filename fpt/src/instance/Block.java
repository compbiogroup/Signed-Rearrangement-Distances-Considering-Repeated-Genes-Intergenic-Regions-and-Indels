package instance;

import java.util.*;
import java.lang.Math;
import main.MCSP;

public class Block {
    ArrayList<Integer> genes;
    ArrayList<Integer> irs;
    boolean inv;
    boolean to_delete;

    public static class InvBlock {
        ArrayList<Integer> genes;
        ArrayList<Integer> irs;
        int sign;
        boolean to_delete;

        private int compareLists(ArrayList<Integer> o1, ArrayList<Integer> o2) {
            for (int i = 0; i < Math.min(o1.size(), o2.size()); i++) {
              int c = o1.get(i).compareTo(o2.get(i));
              if (c != 0) {
                return c;
              }
            }
            return Integer.compare(o1.size(), o2.size());
          }

        public InvBlock(Block b) {
            genes = b.genes;
            irs = b.irs;
            to_delete = b.to_delete;
            if (b.inv) {
                Collections.reverse(genes);
                Collections.reverse(irs);
            }

            ArrayList<Integer> rgenes = new ArrayList<Integer>();
            for (Integer i : genes) {
                rgenes.add(-i);
            }
            Collections.reverse(rgenes);
            ArrayList<Integer> rirs = new ArrayList<Integer>(irs);
            Collections.reverse(rirs);

            int cmp = compareLists(genes,rgenes);
            if (cmp < 0) {
                sign = -1;
            } else if (cmp > 0) {
                sign = 1;
            } else {
                if (compareLists(irs,rirs) <= 0) {
                    sign = 1;
                } else {
                    sign = -1;
                }
            }
            if (sign == -1) {
                genes = rgenes;
                irs = rirs;
            }
        }

        public int size() {
            return genes.size();
        }

        public boolean isNegative() {
            return sign == -1;
        }

        public boolean isDeleted() {
            return to_delete;
        }

        public Integer getIR(int i) {
            return irs.get(i);
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((genes == null) ? 0 : genes.hashCode());
            result = prime * result + ((irs == null) ? 0 : irs.hashCode());
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj)
                return true;
            if (obj == null)
                return false;
            if (getClass() != obj.getClass())
                return false;
            Block.InvBlock other = (Block.InvBlock) obj;
            if (!genes.equals(other.genes))
                return false;
            if (!irs.equals(other.irs))
                return false;
            return true;
        }

        public void print() {
            MCSP.print(genes.toString());
            MCSP.print(irs.toString());
            MCSP.print("sign: " + Integer.toString(sign));
        }

    }

	public Block(boolean to_delete){
        this(to_delete, false);
    }

	public Block(boolean to_delete, boolean inv){
        this.to_delete = to_delete;
        genes = new ArrayList<Integer>();
        irs = new ArrayList<Integer>();
        this.inv = inv;
	}

	public void addGene(Integer gene){
		this.genes.add(gene);
	}

	public void addIR(Integer ir){
        this.irs.add(ir);
    }
    
}
