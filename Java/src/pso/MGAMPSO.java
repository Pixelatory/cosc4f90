package pso;

import base.MGPSO;
import org.kamranzafar.commons.cloner.ObjectCloner;
import util.*;

import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

public class MGAMPSO extends MGPSO {
    private int numOfInfeasibleSols = 0;
    private ArrayList<Triplet<double[], int[][], Double>> gBest;
    private ArrayList<Triplet<double[], int[][], Double>> sArchive;

    public MGAMPSO(String[] seq,
                   int[] n,
                   double w,
                   double c1,
                   double c2,
                   double c3,
                   double vmax,
                   int vmaxiterlimit,
                   double[] term,
                   int maxIter,
                   FitnessFunction[] f,
                   Operator[] ops) throws Exception {
        super(seq, n, w, c1, c2, c3, vmax, vmaxiterlimit, term, maxIter, f, ops);
    }

    public void startPSO() {
        // Begin checks for trivial errors
        for (int amount : n)
            if (amount < 1)
                throw new IllegalArgumentException("Swarm size cannot be < 1");
        if (seq.length < 2)
            throw new IllegalArgumentException("Number of sequences < 2");
        else if (maxIter < 1)
            throw new IllegalArgumentException("Maximum iterations cannot be < 1");

        int highestN = 0;
        int sumN = 0;
        for (int amount : n) {
            if (amount > highestN)
                highestN = amount;

            sumN += amount;
        }

        int colLength = Helper.getColLength(seq);
        int numOfSeqs = seq.length;

        // Initialize main data containers
        double[][][] pPositions = new double[f.length][][];
        double[][][] pPersonalBests = new double[f.length][][];
        double[][] pFitnesses = new double[f.length][];
        int[][][][] pBitStrings = new int[f.length][][][];
        double[][][] pVelocities = new double[f.length][][];
        gBest = new ArrayList<>();
        sArchive = new ArrayList<>();
        double l = 0; // lambda coefficient

        numOfInfeasibleSols = 0;

        // Initialize particles of each sub-swarm
        for (int i = 0; i < f.length; i++) {
            double[][] newPositions = new double[n[i]][4];
            double[][] newVelocities = new double[n[i]][4];
            int[][][] newBitStrings = new int[n[i]][numOfSeqs][colLength];

            for (int j = 0; j < n[i]; j++) {
                double[] tmpPos = new double[4];
                double[] tmpVel = new double[4];
                int[][] bitstring;

                // Velocity is initialized to 0-vector
                tmpVel[0] = 0;
                tmpVel[1] = 0;
                tmpVel[2] = 0;
                tmpVel[3] = 0;

                // Ensure position and bitstring are feasible
                do {
                    tmpPos[0] = ThreadLocalRandom.current().nextDouble(-0.3, 0.3);
                    tmpPos[1] = ThreadLocalRandom.current().nextDouble(0.5, 10);
                    tmpPos[2] = ThreadLocalRandom.current().nextDouble(0.5, 10);
                    tmpPos[3] = ThreadLocalRandom.current().nextDouble(-0.9, -0.5);

                    bitstring = Helper.genBitMatrix(tmpPos, seq, colLength);
                } while(Helper.infeasible(bitstring, seq, ops));

                newPositions[j] = tmpPos;
                newVelocities[j] = tmpVel;
                newBitStrings[j] = bitstring;
                pFitnesses[i][j] = f[i].calculate(bitstring, seq);
            }

            ObjectCloner<double[][]> cloner = new ObjectCloner<>();
            pPositions[i] = newPositions;
            pVelocities[i] = newVelocities;
            pBitStrings[i] = newBitStrings;
            pPersonalBests[i] = cloner.deepClone(newPositions);
        }

        int iter = 0;
        while(iter < maxIter) {
            // If global best is better than termination criteria for respective sub-swarm, then quit
            for (int i = 0; i < f.length; i++) {
                if (gBest.get(i).getThird() < term[i])
                    return;
            }
            iter++;
        }
    }
}
