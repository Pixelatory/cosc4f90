package base;

import PSOs.MGPSO;
import util.FitnessFunction;
import util.Helper;
import util.Operator;
import util.Pair;

import java.util.ArrayList;

/*
    Multi-Guided Binary PSO

    Uses linearly increasing lambda value
    from 0-1 over iterations.
 */
public class MGBPSO extends MGPSO {
    private int numOfInfeasibleSols = 0;
    private ArrayList<ArrayList<Integer>> gBestPos;

    public MGBPSO(ArrayList<String> seq,
                  int n,
                  double w,
                  double c1,
                  double c2,
                  double c3,
                  double vmax,
                  int vmaxiterlimit,
                  double term,
                  int maxIter,
                  Pair<FitnessFunction, Operator>[] f,
                  ArrayList<Operator> ops) {
        super(seq, n, w, c1, c2, c3, vmax, vmaxiterlimit, term, maxIter, f, ops);
    }

    public void startPSO() {
        // Begin checks for trivial errors
        if (n < 1)
            throw new IllegalArgumentException("Swarm size cannot be < 1");
        else if (seq.size() < 2)
            throw new IllegalArgumentException("Number of sequences < 2");
        else if (maxIter < 1)
            throw new IllegalArgumentException("Maximum iterations cannot be < 1");

        int colLength = Helper.getColLength(seq);
        int numOfSeqs = seq.size();

        // Initialize main data containers
        int[][][][] pPositions = new int[f.length][n][numOfSeqs][colLength];
        int[][][][] pPersonalBests = new int[f.length][n][numOfSeqs][colLength];
        double[][][][] pVelocities = new double[f.length][n][numOfSeqs][colLength];

    }
}
