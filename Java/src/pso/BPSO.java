package pso;

import org.kamranzafar.commons.cloner.ObjectCloner;
import util.Helper;
import util.Operator;
import base.PSO;
import util.Sequences;

import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

public class BPSO extends PSO {
    private final double w1;
    private final double w2;
    private int numOfInfeasibleSols = 0;
    private int[][] gBestPos;
    private double gBestFitness;

    private double fitness(int[][] bitmatrix) {
        return Helper.aggregatedFunction(bitmatrix, seq, w1, w2, true, ops);
    }

    public BPSO(String[] seq,
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
        // Begin checks for trivial errors
        if (n < 1)
            throw new IllegalArgumentException("Swarm size cannot be < 1");
        else if (seq.length < 2)
            throw new IllegalArgumentException("Number of sequences < 2");
        else if (maxIter < 1)
            throw new IllegalArgumentException("Maximum iterations cannot be < 1");
        else if (term == null || term.length == 0)
            throw new IllegalArgumentException("Termination criteria is empty");
        else if (term.length > 1)
            System.err.println("Warning: termination array doesn't need to have more than 1 entry for AMPSO.");

        // Vars just to make things more readable
        int colLength = Helper.getColLength(seq);
        int numOfSeq = seq.length;

        // Initialize main data containers
        gBestPos = new int[numOfSeq][colLength];
        int[][][] pPositions = new int[n][numOfSeq][colLength];
        int[][][] pPersonalBests = new int[n][numOfSeq][colLength];
        double[][][] pVelocities = new double[n][numOfSeq][colLength];
        double[] pFitnesses = new double[n]; // so we only calculate fitness once

        numOfInfeasibleSols = 0;

        // Initializing the global best particle to all 0 integer bitmatrix
        for (int i = 0; i < numOfSeq; i++) {
            for (int j = 0; j < colLength; j++) {
                gBestPos[i][j] = 0;
            }
        }

        gBestFitness = fitness(gBestPos); // Initializing the global best fitness

        // Initializing the particles of the swarm (they will always initially be feasible)
        for (int i = 0; i < n; i++) {
            int[][] position = new int[numOfSeq][colLength];
            double[][] velocity = new double[numOfSeq][colLength];

            for (int j = 0; j < numOfSeq; j++) {
                // First the velocity and position is set to all 0
                for (int k = 0; k < colLength; k++) {
                    velocity[j][k] = 0.0;
                    position[j][k] = 0;
                }

                // Now this transforms the position into something feasible
                for (int k = 0; k < (colLength - seq[j].length()); k++) {
                    while (true) {
                        int randNum = ThreadLocalRandom.current().nextInt(0, colLength);
                        if (position[j][randNum] != 1) {
                            position[j][randNum] = 1;
                            break;
                        }
                    }
                }
            }

            ObjectCloner<int[][]> cloner = new ObjectCloner<>();
            pPositions[i] = position;
            pPersonalBests[i] = cloner.deepClone(position);
            pVelocities[i] = velocity;
            pFitnesses[i] = fitness(position);
        }

        // This is where the iterations begin
        int iter = 0; // iteration counter
        while (iter < maxIter && gBestFitness < term[0]) {

            // Checking of fitness is better than global best for updating
            for (int i = 0; i < n; i++) {
                if (pFitnesses[i] > gBestFitness) {
                    ObjectCloner<int[][]> cloner = new ObjectCloner<>();
                    gBestPos = cloner.deepClone(pPersonalBests[i]);
                    gBestFitness = pFitnesses[i];
                }
            }

            // Update each particle's velocity, position, and personal best
            for (int i = 0; i < n; i++) {
                // r1 and r2 are ~ U(0,1)
                double r1 = ThreadLocalRandom.current().nextDouble(0, 1);
                double r2 = ThreadLocalRandom.current().nextDouble(0, 1);

                for (int j = 0; j < numOfSeq; j++) {
                    for (int k = 0; k < colLength; k++) {
                        pVelocities[i][j][k] = w * pVelocities[i][j][k]
                                + r1 * c1 * (pPersonalBests[i][j][k] - pPositions[i][j][k])
                                + r2 * c2 * (gBestPos[j][k] - pPositions[i][j][k]);

                        if (vmaxiterlimit < iter) {
                            if (pVelocities[i][j][k] > vmax)
                                pVelocities[i][j][k] = vmax;
                            else if (pVelocities[i][j][k] < -vmax)
                                pVelocities[i][j][k] = -vmax;
                        }

                        double probability = Helper.Sigmoid(pVelocities[i][j][k]);
                        pPositions[i][j][k] = ThreadLocalRandom.current().nextDouble() < probability ? 1 : 0;
                    }
                }

                double tmpFitness = fitness(pPositions[i]);

                // solution is infeasible, so increment count
                if (tmpFitness == Double.MIN_VALUE)
                    numOfInfeasibleSols++;
                else if (tmpFitness > pFitnesses[i]) {
                    // current fitness is better than personal best's so update it
                    ObjectCloner<int[][]> cloner = new ObjectCloner<>();
                    pFitnesses[i] = tmpFitness;
                    pPersonalBests[i] = cloner.deepClone(pPositions[i]);
                }
            }

            iter++;
        }
    }

    public static void main(String[] args) throws InterruptedException {
        ArrayList<BPSO> bs = new ArrayList<>();

        Operator[] ops = {Operator.lt};

        int n = 30;
        int maxIter = 5000;
        double w1 = 0.5;
        double w2 = 0.5;

        for (int i = 0; i < 30; i++) {
            BPSO b = new BPSO(Sequences.seq1, n, 0.99, 2, 2, 11, 0, new double[]{Double.MAX_VALUE}, maxIter, w1, w2, ops);
            b.start();
            bs.add(b);
        }

        for (int i = 0; i < 30; i++) {
            bs.get(i).join();
            System.out.println("Test Run:");
            System.out.println((((double) bs.get(i).numOfInfeasibleSols / (n * maxIter)) * 100) + "% infeasible");
            System.out.println(Helper.aggregatedFunction(bs.get(i).gBestPos, Sequences.seq1, w1, w2, true, ops));
        }


    }
}