import base.AM;
import util.Helper;
import util.Operator;
import util.Pair;

import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

public class AMPSO extends AM {
    private final double w1;
    private final double w2;
    private int numOfInfeasibleSols = 0;
    private ArrayList<Integer> gBestPos;

    private double fitness(ArrayList<ArrayList<Integer>> bitmatrix) {
        return Helper.aggregatedFunction(bitmatrix, seq, w1, w2, true, ops);
    }

    public AMPSO(ArrayList<String> seq,
                 int n,
                 double w,
                 double c1,
                 double c2,
                 double vmax,
                 int vmaxiterlimit,
                 double term,
                 int maxIter,
                 double w1,
                 double w2,
                 ArrayList<Operator> ops) {
        super(seq, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, ops);
        this.w1 = w1;
        this.w2 = w2;
    }

    public void run() {
        startPSO();
    }

    public void startPSO() {
        if (n < 1)
            throw new IllegalArgumentException("Swarm size cannot be < 1");
        else if (seq.size() < 2)
            throw new IllegalArgumentException("Number of sequences < 2");
        else if (maxIter < 1)
            throw new IllegalArgumentException("Maximum iterations cannot be < 1");
        else if (maxIter == Integer.MAX_VALUE)
            System.err.println("Warning: maximum iterations is set to integer max value");

        // Note: genInterval will always be in sorted where first is less than second (view AM.java)

        // Initialize data containers
        gBestPos = new ArrayList<>();
        double gBestFitness;
        ArrayList<ArrayList<Integer>> gBestBitString = new ArrayList<>();

        ArrayList<ArrayList<Integer>> pPositions = new ArrayList<>();
        ArrayList<ArrayList<Integer>> pPersonalBests = new ArrayList<>();
        ArrayList<ArrayList<Double>> pVelocities = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Integer>>> pBitStrings = new ArrayList<>();
        ArrayList<Double> pFitnesses = new ArrayList<>(); // so we only calculate fitness once

        numOfInfeasibleSols = 0;

        // More vars just to make things more readable
        int colLength = Helper.getColLength(seq);
        int numOfSeq = seq.size();

        while(true) {
            ArrayList<Double> tmpPos = new ArrayList<>();
            tmpPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));
            tmpPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));
            tmpPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));
            tmpPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));



        }
    }
}
