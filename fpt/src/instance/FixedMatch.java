package instance;

import java.util.*;

public class FixedMatch{
	Marker m1;
	Marker m2;
    Integer prev_ir1;
    Integer prev_ir2;
    Integer pos1;
    Integer pos2;
	Added where;
	public enum Added{DR1,DR2,DR3,End,Unspec};

	FixedMatch(Marker m1, Marker m2){
		this.m1 = m1;
		this.m2 = m2;
		this.where = Added.Unspec;
        prev_ir1 = 0;
        prev_ir2 = 0;
	}
	
	FixedMatch(Marker m1, Marker m2,Added where){
		this.m1 = m1;
		this.m2 = m2;
		this.where = where;
        prev_ir1 = 0;
        prev_ir2 = 0;
	}

    @Override
    public String toString() {
        return this.m1.tag() + " -- " + this.m2.tag() + " " + where;
    }
	

    public Integer getM1Pos() {
        return this.pos1;
    }

    public Integer getM2Pos() {
        return this.pos2;
    }

    public void setM1Pos(int pos) {
        this.pos1 = pos;
    }

    public void setM2Pos(int pos) {
        this.pos2 = pos;
    }

    public Integer getPrevIR1() {
        return this.prev_ir1;
    }

    public Integer getPrevIR2() {
        return this.prev_ir2;
    }

}
