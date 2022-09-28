package main;

import algorithm.ReductionRules;
import instance.*;

public class EvaluateDataReduction {

	public static void evaluate(Instance inst) throws Exception {
		System.out.println("n1\tn2\tn1\'\tn2\'\tmatchReduce");
		Instance testInstance = InstanceFactory.copyInstance(inst);
		int removedMatches = testInstance.reduce();
		System.out.println(inst.getGenome1().size()+"\t"+inst.getGenome2().size()+"\t"+
				testInstance.getGenome1().size()+"\t"+testInstance.getGenome2().size()+"\t"+removedMatches);
		//without Rule1
		testWithoutRule(1, inst);
		testWithoutRule(2, inst);
		testWithoutRule(3, inst);
		testWithoutRule(4, inst);
		testInstance.countProperties();
	}	
	
	private static void testWithoutRule(int i,Instance inst) throws Exception {
		System.out.println("Test Without Rule "+i);
		Instance testInstance = InstanceFactory.copyInstance(inst);
		//switch rule off
		ReductionRules.switchRule(false,i);
		int removedMatches = testInstance.reduce();
		System.out.println(inst.getGenome1().size()+"\t"+inst.getGenome2().size()+"\t"+
				testInstance.getGenome1().size()+"\t"+testInstance.getGenome2().size()+"\t"+removedMatches);
		ReductionRules.switchRule(true,i);
		//switch rule back on		
	}
}
