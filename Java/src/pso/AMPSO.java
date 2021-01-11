package pso;

import base.AM;
import util.ArrayCloner;
import util.Helper;
import util.Operator;
import util.Sequences;

import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

public class AMPSO extends AM {
    private final double w1;
    private final double w2;
    private int numOfInfeasibleSols = 0;
    private double[] gBestPos;
    private int[][] gBestBitString;
    private double gBestFitness;

    private double fitness(int[][] bitmatrix) {
        return Helper.aggregatedFunction(bitmatrix, seq, w1, w2, true, ops);
    }

    public AMPSO(String[] seq,
                 int n,
                 double w,
                 double c1,
                 double c2,
                 double vmax,
                 int vmaxiterlimit,
                 double[] term,
                 int maxIter,
                 double w1,
                 double w2,
                 Operator[] ops) {
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
        else if (seq.length < 2)
            throw new IllegalArgumentException("Number of sequences < 2");
        else if (maxIter < 1)
            throw new IllegalArgumentException("Maximum iterations cannot be < 1");
        else if(term == null || term.length == 0)
            throw new IllegalArgumentException("Termination criteria is empty");
        else if(term.length > 1)
            System.err.println("Warning: termination array doesn't need to have more than 1 entry for AMPSO.");

        // Note: genInterval will always be in sorted where first is less than second (view AM.java)

        //  Vars just to make things more readable
        int colLength = Helper.getColLength(seq);
        int numOfSeqs = seq.length;

        // Initialize data containers
        double[][] pPositions = new double[n][4];
        double[][] pPersonalBests = new double[n][4];
        double[][] pVelocities = new double[n][4];
        int[][][] pBitStrings = new int[n][numOfSeqs][colLength];
        double[] pFitnesses = new double[n]; // so we only calculate fitness once

        numOfInfeasibleSols = 0;

        // TODO: Look at papers at how to initialize positions
        // Initializing the global best particle
        do {
            gBestPos = new double[4];
            gBestPos[0] = ThreadLocalRandom.current().nextDouble(-3, 3);
            gBestPos[1] = ThreadLocalRandom.current().nextDouble(-3, 3);
            gBestPos[2] = ThreadLocalRandom.current().nextDouble(-3, 3);
            gBestPos[3] = ThreadLocalRandom.current().nextDouble(-3, 3);

            gBestBitString = Helper.genBitMatrix(gBestPos, seq, colLength);
        } while (Helper.infeasible(gBestBitString, seq, ops));

        gBestFitness = fitness(gBestBitString);

        // Initializing the particles of swarm
        for (int i = 0; i < n; i++) {
            double[] tmpPos = new double[4];
            double[] tmpVel = new double[4];

            tmpPos[0] = ThreadLocalRandom.current().nextDouble(-3, 3);
            tmpPos[1] = ThreadLocalRandom.current().nextDouble(-3, 3);
            tmpPos[2] = ThreadLocalRandom.current().nextDouble(-3, 3);
            tmpPos[3] = ThreadLocalRandom.current().nextDouble(-3, 3);

            tmpVel[0] = 0;
            tmpVel[1] = 0;
            tmpVel[2] = 0;
            tmpVel[3] = 0;


            int[][] bitstring = Helper.genBitMatrix(tmpPos, seq, colLength);

            pPositions[i] = tmpPos;
            pPersonalBests[i] = tmpPos.clone();
            pVelocities[i] = tmpVel;
            pBitStrings[i] = bitstring;
            pFitnesses[i] = fitness(bitstring);
        }

        // This is where the iterations begin
        int iter = 0;
        while (iter < maxIter && gBestFitness < term[0]) {

            // Update global best if possible
            for (int i = 0; i < n; i++) {
                if (Helper.infeasible(pBitStrings[i], seq, ops))
                    numOfInfeasibleSols++;
                else if (pFitnesses[i] > gBestFitness) {
                    gBestPos = pPersonalBests[i].clone();
                    gBestBitString = ArrayCloner.deepcopy(pBitStrings[i]);
                    gBestFitness = pFitnesses[i];
                }
            }

            // Update each particle's velocity, position, and personal best
            for (int i = 0; i < n; i++) {
                double r1 = ThreadLocalRandom.current().nextDouble();
                double r2 = ThreadLocalRandom.current().nextDouble();

                // Velocity update
                for (int j = 0; j < 4; j++) {
                    pVelocities[i][j] = w * pVelocities[i][j] +
                            r1 * c1 * (pPersonalBests[i][j] - pPositions[i][j]) +
                            r2 * c2 * (gBestPos[j] - pPositions[i][j]);

                    if (vmaxiterlimit < iter) {
                        if (pVelocities[i][j] > vmax)
                            pVelocities[i][j] = vmax;
                        else if (pVelocities[i][j] < -vmax)
                            pVelocities[i][j] = -vmax;
                    }

                    // Position update
                    pPositions[i][j] += pVelocities[i][j];
                }

                // Create bitstring from updated position
                int[][] bitstring = Helper.genBitMatrix(pPositions[i], seq, colLength);

                double tmpScore = fitness(bitstring);

                // Update personal best if possible
                if(tmpScore > pFitnesses[i]) {
                    pPersonalBests[i] = pPositions[i].clone();
                    pBitStrings[i] = bitstring;
                    pFitnesses[i] = tmpScore;
                }
            }
            iter++;
        }
    }

    public static void main(String[] args) throws InterruptedException {
        ArrayList<AMPSO> ams = new ArrayList<>();

        Operator[] ops = {Operator.lt};

        int n = 30;
        int maxIter = 5000;
        double w1 = 0.5;
        double w2 = 0.5;

        for (int i = 0; i < 30; i++) {
            AMPSO a = new AMPSO(Sequences.seq1, n, 0.99, 2, 2, 11, 0, new double[]{Double.MAX_VALUE}, maxIter, w1, w2, ops);
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
