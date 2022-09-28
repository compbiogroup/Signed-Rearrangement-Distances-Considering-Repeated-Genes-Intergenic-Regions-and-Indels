package instance;


import java.util.*;

import algorithm.*;
import main.*;

public class Instance {
	Genome genome1,genome2,fullGenome1,fullGenome2;
	ReductionRules reducer;
	// Solver solver;
	TempSol tempSol;
	

	// private int k;
	
	public Instance(){
		this.reducer = new ReductionRules(this);
		// this.solver = new Solver(this);
		this.genome1 = new Genome(1);
		this.genome2 = new Genome(2);
		this.fullGenome1 = new Genome(1);
		this.fullGenome2 = new Genome(2);
		this.tempSol = new TempSol();
	}
	
	public Instance(Genome g1, Genome g2){
		this.reducer = new ReductionRules(this);
		// this.solver = new Solver(this);
		this.genome1 = new Genome(1);
		this.genome2 = new Genome(2);
		this.fullGenome1 = g1;
		this.fullGenome2 = g2;
		this.tempSol = new TempSol();
	}
	
    public void freezeGenomes(){ 
		for (Iterator<Marker> it = genome1.iterator(); it.hasNext();) {
			Marker current = it.next();
			fullGenome1.addMarker(new Marker(current,genome1));
		}
		for (Iterator<Marker> it = genome2.iterator(); it.hasNext();) {
			Marker current = it.next();
			fullGenome2.addMarker(new Marker(current,genome2));			
		}
    }
	
	public int reduce() throws Exception {
		return reducer.applyRules();
	}
	
	public void print(){
		// System.out.println("Remaining k: "+k);
		genome1.info();
		genome2.info();
	}
	
	public void printFull(){
		// System.out.println("Remaining k: "+k);
		genome1.print();
		genome2.print();
		// fullGenome1.print();
		// fullGenome2.print();
	}
	
	// public void foundBlock(){
	// 	this.k--;
	// }
	

	
	
	public Genome getGenome1(){
		return genome1;
	}
	
	public Genome getGenome2(){
		return genome2;
	}

	public Genome getFullGenome1(){
		return fullGenome1;
	}
	
	public Genome getFullGenome2(){
		return fullGenome2;
	}

	//check whether a Marker is parallel to some other fixed marker
	public boolean checkParallel(Marker m){		
		MCSP.print("Checking whether new Fixed Edge becomes parallel");	
		for (Marker m2: genome1.markers){
			if (checkParallel(m,m2)) return true;
		}		
		return false;
	}
	//check whether two fixed Genes in Genome 1 are parallel, first Argument must be before second Argument
	public boolean checkParallel(Marker m1, Marker m2){		
		//basic checks
		if (!(genome1.contains(m1)&&genome1.contains(m2))) return false;
		if (!(m1.isFixed()&&m2.isFixed())) return false;
		if (m1.getChromosome()!=m2.getChromosome()) return false;		
		MCSP.printV("Checking parallel for fixed markers on same genome");		
		//check order of markers
		Marker first,second;
		if (genome1.getMarkerIndex(m1)<genome1.getMarkerIndex(m2)) {
			first = m1; second = m2;
		}
		else if (genome1.getMarkerIndex(m1)>genome1.getMarkerIndex(m2)) {
			first = m2; second = m1;
		}		
		else {
			MCSP.printV("Both are the same marker");
			return false; //both are the same marker
		}
		Marker firstMatch = first.getMatch();
		Marker secondMatch = second.getMatch();
		if (firstMatch.getChromosome()!=secondMatch.getChromosome()) return false;
		MCSP.printV("Matches on Same Chromosome");
		if (first.isForwardMatch(firstMatch,true)&&second.isForwardMatch(secondMatch,false)){
			MCSP.printV("Both forward");
			//check forward
			if (genome2.getMarkerIndex(secondMatch)<=genome2.getMarkerIndex(firstMatch)) return false;
			MCSP.printV("m2 after m1");
			ListIterator<Marker> gIt = genome1.iterator(genome1.getMarkerIndex(first));
			ListIterator<Marker> mIt = genome2.iterator(genome2.getMarkerIndex(firstMatch));
			Marker gPoint= gIt.next();
			Marker mPoint= mIt.next();			
			while (gIt.hasNext()&&mIt.hasNext()){
				gPoint= gIt.next();
				mPoint= mIt.next();
				if (gPoint.equals(second)&&mPoint.equals(secondMatch)){
					MCSP.printV("Success");
					return true;
				}
				if (!gPoint.isForwardMatch(mPoint,false)){
					MCSP.printV("No Success");
					return false;
				}
			}
			MCSP.print("Should not be reached");			
		}
		else if (first.isBackwardMatch(firstMatch,true)&&second.isBackwardMatch(secondMatch,false)){
			//check backward
			if (genome2.getMarkerIndex(secondMatch)>=genome2.getMarkerIndex(firstMatch)) return false;
			MCSP.printV("m2 before m1");
			ListIterator<Marker> gIt = genome1.iterator(genome1.getMarkerIndex(first));
			ListIterator<Marker> mIt = genome2.iterator(genome2.getMarkerIndex(firstMatch));
			Marker gPoint= gIt.next();
			Marker mPoint= mIt.next();			
			mPoint= mIt.previous();			
			while (gIt.hasNext()&&mIt.hasPrevious()){
				gPoint= gIt.next();
				mPoint= mIt.previous();
				if (gPoint.equals(second)&&mPoint.equals(secondMatch)){
					MCSP.printV("Success");
					return true;
				}
				if (!gPoint.isBackwardMatch(mPoint,false)){
					MCSP.printV("No Success");
					return false;
				}
				MCSP.printV("One match");
			}
			MCSP.print("Should not be reached");
		}
		return false;
	}


	// public int getK() {
	// 	return k;
	// }


	// public void setK(int k) {
	// 	this.k = k;
	// }


	public void makeMatchesConsistent() {
		cleanUpGenome(genome1);
		cleanUpGenome(genome2);
	}

	public Marker getFromGenome1(Marker m, Marker m2) {
		if (genome1.contains(m)) return m;
		if (genome1.contains(m2)) return m2;
		else return null;
	}
	
	private void cleanUpGenome(Genome g){
		for (Iterator<Marker> it = g.iterator(); it.hasNext();) {
			Marker m = it.next();
			for (Iterator<Marker> itmatches = m.getMatches().iterator(); itmatches.hasNext();) {
				Marker m2 = itmatches.next();
				if (!m2.getMatches().contains(m)){
					itmatches.remove();
				}
			}
		}
	}
	
	public void removeNegative(){
		MCSP.print("Removing Negative");
		genome1.removeNegative();
		genome2.removeNegative();
		// genome1.removeNonMatched();
		// genome2.removeNonMatched();
		
	}
	
	public Genome getOther(Genome g){
		if (genome1==g) return genome2;
		else return genome1;
	}

	public TempSol getTempSol() {
		return tempSol;
	}

	public void setTempSol(TempSol tempSol) {
		this.tempSol = tempSol;
	}
	
	public void countProperties(){
		int deg2=0;
		int clust2=0;
		for (Marker m: genome1.markers){
			Marker succ = genome1.getSuccessor(m);			
			if (m.matchNum()==3){
				deg2++;

			}
			if (m.matchNum()>=1&&succ!=null&&succ.isMatch(m.getMatch())) clust2++;
		}
		System.out.println("Degree three in genome1 "+deg2);
		System.out.println("Two cluster "+clust2);
		
		deg2=0;
		clust2=0;
		for (Marker m: genome2.markers){
			Marker succ = genome2.getSuccessor(m);
			if (m.matchNum()==3){
				deg2++;
			}
			if (m.matchNum()>=1&&succ!=null&&succ.isMatch(m.getMatch())) clust2++;

		}
		System.out.println("Degree three in genome2 "+deg2);
		System.out.println("Two cluster "+clust2);
	}
	
}
