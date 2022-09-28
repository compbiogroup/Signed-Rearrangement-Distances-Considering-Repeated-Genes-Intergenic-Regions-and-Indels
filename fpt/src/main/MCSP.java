package main;

import instance.*;
import algorithm.Solver;
import utilities.Timer;
import java.util.ArrayList;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.BufferedWriter;

public class MCSP {
    public static final boolean verbose = false;
    public static final boolean veryVerbose = false;
    public static boolean commonCost = false;
    public static boolean parallel = true;
    public static boolean reduce = true;


    /**
     * @param args
     */
    public static void main(String[] args) throws Exception {

        ArrayList<String[]> line_groups = new ArrayList<String[]>();
        System.out.println("Input File: " + args[0]);
        if (args.length < 3) {
            throw new IllegalArgumentException("Invalid number of command line arguments.");
        }
        if (args[2].equals("CC")) {
            commonCost = true;
        } else if (args[2].equals("NC")) {
            commonCost = false;
        } else {
            throw new IllegalArgumentException("Invalid command line argument (must be CC or NC).");
        }

        // ArrayList<Integer> missing = new ArrayList<Integer>();
        // try {
        //     BufferedReader inFile = new BufferedReader(new FileReader(args[3]));
        //     while (true) {
        //         String line = inFile.readLine();
        //         if (line == null)
        //             break;
        //         missing.add(Integer.parseInt(line));
        //     }
        //     inFile.close();
        // } catch (Exception e) {
        //     System.out.println("Exception in main: " + e.getMessage());
        // }

        try {
            BufferedReader inFile = new BufferedReader(new FileReader(args[0]));
            int id = 1;
            while (true) {
                String[] lines = new String[5];
                lines[0] = inFile.readLine();
                lines[1] = inFile.readLine();
                lines[2] = inFile.readLine();
                lines[3] = inFile.readLine();
                lines[4] = args[1] + "_" + String.format("%04d" , id);
                if (lines[3] == null)
                    break;
                // if (missing.contains(id)) {
                line_groups.add(lines);
                // }
                id++;
            }
            inFile.close();
        } catch (Exception e) {
            System.out.println("Exception in main: " + e.getMessage());
        }

        if (parallel) {
            line_groups.parallelStream().forEach(x -> runOneInstance(x));
        } else {
            line_groups.stream().forEach(x -> runOneInstance(x));
        }
    }

    public static void runOneInstance(String[] lines) {
        Instance inst;
        inst = InstanceFactory.createInstancefromLines(lines[0], lines[1], lines[2], lines[3]);
        Solver solver = new Solver();
        long startTime = System.currentTimeMillis();
        Timer timer = new Timer();
        Boolean finish = true;
        try {
            TempSol.ReducedStr sol = null;
            try {
                sol = solver.solve(inst, timer);
            } catch (InterruptedException e) {
                sol = solver.getSol();
                finish = false;
            }
            if (sol == null) {
                throw new Exception("No solution found.");
            } else {
                BufferedWriter outFile = new BufferedWriter(new FileWriter(lines[4]));
                outFile.write(sol.str);
                outFile.newLine();
                outFile.write("# Time (s): " + timer.elapsed_time());
                outFile.newLine();
                if (!finish) {
                    outFile.write("# Not Exact");
                    outFile.newLine();
                }
                outFile.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Exception in runOneInstanc: " + e.getMessage());
        }
    }

    public static void print(String s) {
        if (verbose)
            System.out.println(s);
    }

    public static void printV(String s) {
        if (veryVerbose)
            System.out.println(s);
    }

    public static void checkCorrectness(boolean quitIfWrong, TempSol sol, Instance inst) {
        try {
            sol.checkIsComplete(inst);
            MCSP.print("Solution checked: correct.");
        } catch (Exception e) {
            MCSP.print("Solution checking returned an error!");
            MCSP.print("[Instance : ]");
            inst.printFull();
            MCSP.print("[Proposed solution : ]");
            sol.print();
            if (quitIfWrong)
                throw new Error(e);
            else
                e.printStackTrace();

        }
    }

}
