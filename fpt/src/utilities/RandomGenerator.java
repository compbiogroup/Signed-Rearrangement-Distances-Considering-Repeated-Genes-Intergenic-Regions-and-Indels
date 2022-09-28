package utilities;
import java.util.ArrayList;
import java.util.Arrays;


public class RandomGenerator {

	/**
	 * @param args
	 */
	public int [] genome0;
	public int [] genome1;
	public int [] sgn0;
	public int [] sgn1;
	public int [] ir0;
	public int [] ir1;
	
	int numOfGeneFamilies;
	int maxD;
	int n;
	int k;
	int totalNoise;
	
	

	public RandomGenerator(int n, int numOfGeneFamilies, int maxD, 
			int originalK, int totalNoise) throws Exception {
		super();
		this.numOfGeneFamilies = numOfGeneFamilies;
		this.maxD = maxD;
		this.n = n;
		this.k = originalK;
		this.totalNoise = totalNoise;
	


	// working part: generates the two genomes given the complete list of parameters.	


		int i;
		
		genome0 = new int[n];
		genome1 = new int[n];
		sgn0 = new int[n];
		sgn1 = new int[n];
		ir0 = new int[n];
		ir1 = new int[n];
		
		//check parameters
		if (n+totalNoise> maxD*numOfGeneFamilies) throw new Exception("generator error: n and totalNoise are too big compared to d and number of families.");		
		if (n< totalNoise*3) throw new Exception("generator error: n is too small compared to noise.");		
		if (n<1||k<1) throw new Exception("n and k must be strictly positive");
		
		
		//fill signs (all positive here !)
		for (i=0; i<n; i++) {
			sgn0[i]=+1;		
			sgn1[i]=+1;
		}

		//fill intergenic regions
		for (i=0; i<n-1; i++) {
			ir0[i]=0;
			ir1[i]=0;
		}
		
		//create maxD available copies of each family
		ArrayList<Integer> availableCopies= new ArrayList<Integer>();
		for (i=0; i<numOfGeneFamilies; i++) {
			availableCopies.add(maxD);			
		}

		//initialize noisy parts
		ArrayList<ArrayList<Integer>> noise0 = new ArrayList<ArrayList<Integer>> ();
		ArrayList<ArrayList<Integer>> noise1 = new ArrayList<ArrayList<Integer>> ();
		for (i=0; i<=k; i++) { // k+1 iterations
			noise0.add(new ArrayList<Integer>());
			noise1.add(new ArrayList<Integer>());
		}
		
		//fill noisy parts
		
		for (i=0; i<totalNoise; i++) {
			int destinationPart = rand(k+1);
			int geneFamily=getCopy(availableCopies,0);
			noise0.get(destinationPart).add(geneFamily);
			 geneFamily=getCopy(availableCopies,1);
			noise1.get(destinationPart).add(geneFamily);
		}
		
		//initialize blocks
		ArrayList<ArrayList<Integer>> blocks = new ArrayList<ArrayList<Integer>> ();
		for (i=0; i<k; i++) {
			blocks.add(new ArrayList<Integer>());
		}
		
		//fill blocks		
		for (i=0; i<n-totalNoise; i++) {
			int destinationBlock = rand(k);
			int geneFamily=getCopy(availableCopies,2);
			blocks.get(destinationBlock).add(geneFamily);			
		}
		
		//permutation for string 2
		ArrayList<Integer> freeBlocks= new ArrayList<Integer>();
		for (i=0; i<k; i++) {
		   freeBlocks.add(i);
		}
		
		// create strings
		
		int pos0=addArray(genome0, noise0.get(0),0);	
		int pos1=addArray(genome1, noise1.get(0),0);	
		for (i=k-1; i>=0; i--) { //decreasing i!
			pos0=addArray(genome0, blocks.get(i),pos0);
			int j= rand(i+1);			
			pos1=addArray(genome1, blocks.get(freeBlocks.remove(j)),pos1);
			pos0=addArray(genome0, noise0.get(i+1),pos0);	
			pos1=addArray(genome1, noise1.get(i+1),pos1);
		}		
				
	}
	

	//random utility
	private int rand(int x)
	{
		return (int) Math.floor(Math.random()*x);
	}
	
	//  finds a random family which has an available copy, and decrease accordingly the number
	//  of remaining available copies.	
	//  Warning: infinite loop if there are no more available copies.
	//  parity : 0/1 for even/odd only; 2 for all
	private int getCopy(ArrayList<Integer> availableCopies, int parity) {
		int x;		
	    do {
	    	x=rand(numOfGeneFamilies);
	    } while (availableCopies.get(x)==0 || !(parity==2 || x%2 == parity));
	    availableCopies.set(x, availableCopies.get(x)-1);
		return x;
	}

	//add elements for 'genome' taken from 'array'; get and return position in the genome
	
	private int addArray(int [] genome, ArrayList<Integer> array, int position) {

		for (Integer x: array) {
			genome[position]=x;
			position++;			
		}
		return position;
	}
	
	


	public static void main(String[] args) {
		try {
			RandomGenerator G =new RandomGenerator(20, 4, 6, 3, 2);
			System.out.println(Arrays.toString(G.genome0));
			System.out.println(Arrays.toString(G.genome1));
            
		} catch (Exception e) {//should not happen (check parameters!)			
			e.printStackTrace();
		}
	}

}
