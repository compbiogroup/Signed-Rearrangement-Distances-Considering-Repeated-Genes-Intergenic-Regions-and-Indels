package algorithm;

import java.util.*;

import main.MCSP;
import utilities.Timer;
import instance.*;

public class Solver {
    Instance inst = null;
    SampleGraph sg = null;
    TempSol.ReducedStr best_sol = null;
    Timer timer;
    Boolean use_timer;
    Boolean run_del;

    public Solver() {
    }

    public TempSol.ReducedStr solve(Instance inst, Timer timer) throws Exception {
        this.timer = timer;
        this.use_timer = true;
        this.run_del = false;
        if (MCSP.commonCost) {
            rec_solve(inst);
        } else {
            rec_solve(inst);
            this.run_del = true;
            rec_solve(inst);
        }
        return best_sol;
    }

    public TempSol.ReducedStr getSol() {
        return best_sol;
    }

    public void rec_solve(Instance inst) throws Exception {
        if (use_timer && timer.done()) {
            throw new InterruptedException();
        }

        inst.reduce();
        this.inst = inst;
        if (best_sol != null && inst.getTempSol().getK() >= best_sol.cost) {
            return;
        }
        sg = new SampleGraph(inst, best_sol);
        // inst.getGenome1().print();
        // inst.getGenome2().print();
        // sg.print();
        if (sg.isNoInstance) {
            MCSP.print("Too many black edges");
        } else if (sg.hasBlackParallel()) {
            MCSP.print("No because of parallel black edges");
        } else if (best_sol != null
                && inst.getTempSol().getK() + sg.blackEdges.size() - (MCSP.commonCost ? 0 : 1) >= best_sol.cost) {
            MCSP.print("No because too many black edges");
        } else {
            Marker branchMarker = sg.findIsolated();
            if (branchMarker != null) {
                MCSP.print("Found isolated marker " + branchMarker.tag());
                // branch on this marker
                branch(branchMarker);
                return;
            }
            // check for rare odd paths
            Set<Marker> branchMarkers = sg.findOddPath(true);
            if (branchMarkers != null) {
                branch(branchMarkers);
                return;
            }
            if (run_del) {
                branch_del();
            } else {
                noOddPaths();
            }
        }
    }

    private void noOddPaths() throws Exception {
        MCSP.print("Only Cycles and Even Paths");
        Instance local_inst = this.finishSolution();
        TempSol.ReducedStr ans = local_inst.getTempSol().getReducedGenomes(inst.getFullGenome1(),
            inst.getFullGenome2());
        if (ans != null && (best_sol == null || ans.cost < best_sol.cost)) {
            best_sol = ans;
        }
    }

    private void branch_del() throws Exception {
        if (use_timer && timer.done()) {
            throw new InterruptedException();
        }

        Set<Marker> branchMarkers = sg.findOddPath(false);
        if (branchMarkers == null) {
            noOddPaths();
        } else {
            for (Marker m : branchMarkers) {
                if (use_timer && timer.done()) {
                    throw new InterruptedException();
                }
                if (!m.isRare()) {
                    SampleGraph.EdgeData ed = sg.new EdgeData();
                    if (sg.breakOddPath(m, ed)) {
                        branch_del();
                    }
                    sg.restoreOddPath(ed);
                }
            }

        }
    }

    private void branch(Set<Marker> branchMarkers) throws Exception {
        MCSP.print("Branching on Marker Set");
        for (Marker m : branchMarkers) {
            MCSP.print("Branching on Marker:" + m.tag());
            branch(m);
        }
    }

    private void branch(Marker branchMarker) throws Exception {
        Set<Marker> oldMatches = branchMarker.getMatches();
        for (Marker matchMarker : oldMatches) {
            if (use_timer && timer.done()) {
                throw new InterruptedException();
            }
            MCSP.print("Fixing edge to " + matchMarker.tag());
            Set<Marker> matchMarkerMatches = matchMarker.getMatches();

            // make Markers fixed in branches
            Set<Marker> bNewMatch = new HashSet<Marker>();
            bNewMatch.add(matchMarker);
            Set<Marker> mNewMatch = new HashSet<Marker>();
            mNewMatch.add(branchMarker);
            branchMarker.setMatches(bNewMatch);
            matchMarker.setMatches(mNewMatch);
            Marker g1Marker = inst.getFromGenome1(branchMarker, matchMarker);
            if (!inst.checkParallel(g1Marker)) {
                // copy instance
                Instance newInst = InstanceFactory.copyInstance(inst);
                newInst.makeMatchesConsistent();
                Instance tmp = inst;
                rec_solve(newInst);
                inst = tmp;
            } else {
                MCSP.print("Aborted branch (Would give parallel black)");
            }
            // restore matchMarkers Matches
            matchMarker.setMatches(matchMarkerMatches);
        }
        // restore branchMarkers Matches
        branchMarker.setMatches(oldMatches);
    }

    // add all black matches as fixes, all green edges of cycles and green even
    // paths,
    // and all red edges of red even paths
    private Instance finishSolution() throws Exception {
        SampleGraph local_sg = new SampleGraph(sg);
        local_sg.addFixes(timer, use_timer);
        return local_sg.inst;
    }

}
