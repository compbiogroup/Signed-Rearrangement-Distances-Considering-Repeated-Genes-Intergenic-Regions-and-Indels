package utilities;

import java.util.*;

public class Tools {

	public static int[] randomPosArray(int length,int max){
		int[] result = new int[length];
		Random random = new Random();
		for (int i=0;i<length;i++){
			result[i]= random.nextInt(max);
		}
		return result;
	}
	
	public static int[] randomPlusMinusArray(int length){
		int[] result = new int[length];
		Random random = new Random();
		for (int i=0;i<length;i++){
			if(random.nextBoolean()){
				result[i]= 1;
			}
			else {
				result[i]=-1;
			}
		}
		return result;
	}
	
	public static int[] plusArray(int length){
		int[] result = new int[length];
		for (int i=0;i<length;i++){
				result[i]= 1;
		}
		return result;
	}
	
}
