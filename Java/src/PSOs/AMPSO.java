package PSOs;

import base.AM;
import util.Helper;
import util.Operator;
import util.Sequences;

import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

public class AMPSO extends AM {
    private final double w1;
    private final double w2;
    private int numOfInfeasibleSols = 0;
    private ArrayList<Double> gBestPos;
    private ArrayList<ArrayList<Integer>> gBestBitString;

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
        double gBestFitness;

        ArrayList<ArrayList<Double>> pPositions = new ArrayList<>();
        ArrayList<ArrayList<Double>> pPersonalBests = new ArrayList<>();
        ArrayList<ArrayList<Double>> pVelocities = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Integer>>> pBitStrings = new ArrayList<>();
        ArrayList<Double> pFitnesses = new ArrayList<>(); // so we only calculate fitness once

        numOfInfeasibleSols = 0;

        // More vars just to make things more readable
        int colLength = Helper.getColLength(seq);

        // TODO: Look at papers at how to initialize positions
        // Initializing the global best particle
        do {
            gBestPos = new ArrayList<>();
            gBestPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));
            gBestPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));
            gBestPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));
            gBestPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));

            gBestBitString = Helper.genBitMatrix(gBestPos, seq, colLength);
        } while (Helper.infeasible(gBestBitString, seq, ops));

        gBestFitness = fitness(gBestBitString);

        // Initializing the particles of swarm
        for (int i = 0; i < n; i++) {
            ArrayList<Double> tmpPos = new ArrayList<>();
            ArrayList<Double> tmpVel = new ArrayList<>();

            tmpVel.add(0.0);
            tmpVel.add(0.0);
            tmpVel.add(0.0);
            tmpVel.add(0.0);

            tmpPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));
            tmpPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));
            tmpPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));
            tmpPos.add(ThreadLocalRandom.current().nextDouble(-3, 3));

            ArrayList<ArrayList<Integer>> bitstring = Helper.genBitMatrix(tmpPos, seq, colLength);

            pPositions.add(tmpPos);
            pPersonalBests.add(tmpPos);
            pVelocities.add(tmpVel);
            pBitStrings.add(bitstring);
            pFitnesses.add(fitness(bitstring));
        }

        // This is where the iterations begin
        int iter = 0;
        while (iter < maxIter && gBestFitness < term) {

            // Update global best if possible
            for (int i = 0; i < n; i++) {
                if (Helper.infeasible(pBitStrings.get(i), seq, ops)) {
                    numOfInfeasibleSols++;
                } else if (pFitnesses.get(i) > gBestFitness) {
                    gBestPos = (ArrayList<Double>) pPositions.get(i).clone();
                    gBestBitString = Helper.copyArray(pBitStrings.get(i));
                    gBestFitness = pFitnesses.get(i);
                }
            }

            // Update each particle's velocity, position, and personal best
            for (int i = 0; i < n; i++) {
                double r1 = ThreadLocalRandom.current().nextDouble();
                double r2 = ThreadLocalRandom.current().nextDouble();

                // Velocity update
                for (int j = 0; j < 4; j++) {
                    pVelocities.get(i).set(j, w * pVelocities.get(i).get(j) +
                            r1 * c1 * (pPersonalBests.get(i).get(j) - pPositions.get(i).get(j)) +
                            r2 * c2 * (gBestPos.get(j) - pPositions.get(i).get(j)));

                    if (vmaxiterlimit < iter) {
                        if (pVelocities.get(i).get(j) > vmax)
                            pVelocities.get(i).set(j, vmax);
                        else if (pVelocities.get(i).get(j) < -vmax) {
                            pVelocities.get(i).set(j, -vmax);
                        }
                    }

                    // Position update
                    pPositions.get(i).set(j, pPositions.get(i).get(j) + pVelocities.get(i).get(j));
                }

                // Create bitstring from updated position
                ArrayList<ArrayList<Integer>> bitstring = Helper.genBitMatrix(pPositions.get(i), seq, colLength);

                double tmpScore = fitness(bitstring);

                // Update personal best if possible
                if(tmpScore > pFitnesses.get(i)) {
                    pPersonalBests.set(i, (ArrayList<Double>) pPositions.get(i).clone());
                    pBitStrings.set(i, bitstring);
                    pFitnesses.set(i, tmpScore);
                }
            }
            iter++;
        }
    }

    public static void main(String[] args) throws InterruptedException {
        ArrayList<AMPSO> ams = new ArrayList<>();

        ArrayList<Operator> ops = new ArrayList<>();
        ops.add(Helper.lt);

        int n = 30;
        int maxIter = 5000;
        double w1 = 0.5;
        double w2 = 0.5;

        for (int i = 0; i < 30; i++) {
            AMPSO a = new AMPSO(Sequences.seq1, n, 0.99, 2, 2, 11, 0, Double.MAX_VALUE, maxIter, w1, w2, ops);
            a.start();
            ams.add(a);
        }

        for (int i = 0; i < 30; i++) {
            ams.get(i).join();
            System.out.println("Test Run:");
            System.out.println((((double) ams.get(i).numOfInfeasibleSols / (n * maxIter)) * 100) + "% infeasible");
            System.out.println(Helper.aggregatedFunction(ams.get(i).gBestBitString, Sequences.seq1, w1, w2, true, ops));
        }
    }
}
