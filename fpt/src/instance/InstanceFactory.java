package instance;

import utilities.Tools;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Arrays;
import java.lang.Math;

import utilities.RandomGenerator;

import main.MCSP;

public class InstanceFactory {

	public static Instance createInstanceFromArrays(int[] g1, int[] g2, int[] s1, int[] s2, int[] i1, int[] i2){
		Instance inst = new Instance();
		Genome genome1 = inst.getGenome1();		
		Genome genome2 = inst.getGenome2();
        int id = 0;

        genome1.addMarker(new Marker(id, s1[0], 0, g1[0],genome1, 0, i1[0]));			
        id++;
		for (int i = 1; i < g1.length - 1; i++) {
			genome1.addMarker(new Marker(id, s1[i], i, g1[i],genome1, i1[i-1], i1[i]));			
            id++;
		}
        genome1.addMarker(new Marker(id, s1[g1.length - 1], g1.length - 1, g1[g1.length - 1],genome1, i1[g1.length - 2], 0));			
        id++;
        genome2.addMarker(new Marker(id, s2[0], 0, g2[0],genome2, 0, i2[0]));			
        id++;
		for (int i = 1; i < g2.length - 1; i++) {
			genome2.addMarker(new Marker(id, s2[i], i, g2[i],genome2, i2[i-1], i2[i]));			
            id++;
		}
        genome2.addMarker(new Marker(id, s2[g2.length - 1], g2.length - 1, g2[g2.length - 1],genome2, i2[g2.length - 2], 0));			
		for (Iterator<Marker> it1 = genome1.iterator(); it1.hasNext();) {
			Marker m1 = it1.next();
			for (Iterator<Marker> it2 = genome2.iterator(); it2.hasNext();) {
				Marker m2 = it2.next();
					if (m1.name.equals(m2.name)){
						m1.addMatch(m2);
						m2.addMatch(m1);
					}
			}
		}
		genome1.testGenome(genome2);
		genome2.testGenome(genome1);
        inst.freezeGenomes();
		createGenomePointers(inst.getFullGenome1(), inst.getFullGenome2());
		createGenomePointers(inst.getFullGenome2(), inst.getFullGenome1());
		return inst;
	}
	
	public static Instance createInstancefromLines(String l1, String l2, String l3, String l4){
        int[] g1, g2,i1, i2, s1, s2;
        l1 = Integer.toString(0) + " " + l1 + " " + Integer.toString(10000000);
        l3 = Integer.toString(0) + " " + l3 + " " + Integer.toString(10000000);
        g1 = (Arrays.asList(l1.split(" "))).stream().mapToInt(num -> Math.abs(Integer.parseInt(num))).toArray();
        i1 = (Arrays.asList(l2.split(" "))).stream().mapToInt(Integer::parseInt).toArray();
        s1 = (Arrays.asList(l1.split(" "))).stream().mapToInt(num -> Integer.signum(Integer.parseInt(num))).toArray();
        g2 = (Arrays.asList(l3.split(" "))).stream().mapToInt(num -> Math.abs(Integer.parseInt(num))).toArray();
        s2 = (Arrays.asList(l3.split(" "))).stream().mapToInt(num -> Integer.signum(Integer.parseInt(num))).toArray();
        i2 = (Arrays.asList(l4.split(" "))).stream().mapToInt(Integer::parseInt).toArray();
		return createInstanceFromArrays(g1, g2, s1, s2, i1, i2);
	}
	
	public static void createGenomePointers(Genome g1, Genome g2){
		if (MCSP.verbose) System.out.println("Creating Gene Pointers");
		for (Iterator<Marker> it = g1.geneList().iterator(); it.hasNext();){
			Marker gene1 = it.next();
			for (Iterator<Integer> itid = gene1.matchIds.iterator(); itid.hasNext();){
				gene1.addMatch(g2.getMarkerById(itid.next()));
			}
		}
	}
	
	public static Instance copyInstance(Instance inst){
		Instance result = new Instance(inst.getFullGenome1(), inst.getFullGenome2());
		Genome newGenome1 = result.getGenome1();
		Genome newGenome2 = result.getGenome2();		
		Genome oldGenome1 = inst.getGenome1();
		Genome oldGenome2 = inst.getGenome2();
		for (Iterator<Marker> it = oldGenome1.iterator(); it.hasNext();) {
			Marker current = it.next();
			newGenome1.addMarker(new Marker(current,newGenome1));
			
		}
		for (Iterator<Marker> it = oldGenome2.iterator(); it.hasNext();) {
			Marker current = it.next();
			newGenome2.addMarker(new Marker(current,newGenome2));			
		}		
		createGenomePointers(newGenome1, newGenome2);
		createGenomePointers(newGenome2, newGenome1);
		result.setTempSol(inst.getTempSol().copyTempSol());
		return result;
	}
	
}
