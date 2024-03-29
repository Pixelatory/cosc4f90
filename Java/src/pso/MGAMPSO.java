package pso;

import base.MGPSO;
import org.kamranzafar.commons.cloner.ObjectCloner;
import util.*;

import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
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

    public void run() {
        try {
            startPSO();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void startPSO() throws Exception {
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
        double[][][] pPositions = new double[f.length][highestN][4];
        double[][][] pPersonalBests = new double[f.length][highestN][4];
        double[][][] pVelocities = new double[f.length][highestN][4];
        double[][] pFitnesses = new double[f.length][highestN];
        int[][][][] pBitStrings = new int[f.length][highestN][numOfSeqs][colLength];
        gBest = new ArrayList<>();
        sArchive = new ArrayList<>();
        double l = 0; // lambda coefficient

        // Initialize cloners
        ObjectCloner<double[]> positionCloner = new ObjectCloner<>();
        ObjectCloner<int[][]> bitmatrixCloner = new ObjectCloner<>();

        numOfInfeasibleSols = 0;

        // Initialize particles of each sub-swarm
        for (int i = 0; i < f.length; i++) {

            for (int j = 0; j < n[i]; j++) {

                // Velocity is initialized to 0-vector
                pVelocities[i][j][0] = 0;
                pVelocities[i][j][1] = 0;
                pVelocities[i][j][2] = 0;
                pVelocities[i][j][3] = 0;

                int[][] bitstring;

                do {
                    // Initialize position to random double in (-1, 1) for all positions
                    pPositions[i][j][0] = ThreadLocalRandom.current().nextDouble(-1, 1);
                    pPositions[i][j][1] = ThreadLocalRandom.current().nextDouble(-1, 1);
                    pPositions[i][j][2] = ThreadLocalRandom.current().nextDouble(-1, 1);
                    pPositions[i][j][3] = ThreadLocalRandom.current().nextDouble(-1, 1);

                    bitstring = Helper.genBitMatrix(pPositions[i][j], seq, colLength);
                } while (Helper.infeasible(bitstring, seq, ops)
                        || pPositions[i][j][1] * pPositions[i][j][2] == 0);

                pPersonalBests[i][j] = positionCloner.deepClone(pPositions[i][j]);
                pBitStrings[i][j] = bitstring;
                pFitnesses[i][j] = f[i].calculate(bitstring, seq);
            }
            gBest.add(new Triplet<>(positionCloner.deepClone(pPositions[i][0]),
                    bitmatrixCloner.deepClone(pBitStrings[i][0]),
                    pFitnesses[i][0]));
        }

        int iter = 0;
        while (iter < maxIter) {
            // If global best is better than termination criteria for respective sub-swarm, then quit
            for (int i = 0; i < f.length; i++) {
                if (gBest.get(i).getThird() < term[i])
                    return;
            }

            // Update personal best, global best and archive if applicable
            for (int i = 0; i < f.length; i++) {
                for (int j = 0; j < n[i]; j++) {
                    if (Helper.infeasible(pBitStrings[i][j], seq, ops)
                            || (pPositions[i][j][1] * pPositions[i][j][2]) == 0)
                        numOfInfeasibleSols++;
                    else {
                        // Add to archive if possible
                        Helper.addToArchiveA(seq, sArchive, new Pair<>(pPositions[i][j], pBitStrings[i][j]), f, sumN);

                        // Update global best if possible
                        if (pFitnesses[i][j] < gBest.get(i).getThird()) {
                            gBest.get(i).setFirst(positionCloner.deepClone(pPersonalBests[i][j]));
                            gBest.get(i).setSecond(bitmatrixCloner.deepClone(pBitStrings[i][j]));
                            gBest.get(i).setThird(pFitnesses[i][j]);
                        }
                    }
                }
            }

            // Update velocity and position of particles
            for (int i = 0; i < f.length; i++) {
                for (int j = 0; j < n[i]; j++) {
                    double r1 = ThreadLocalRandom.current().nextDouble();
                    double r2 = ThreadLocalRandom.current().nextDouble();
                    double r3 = ThreadLocalRandom.current().nextDouble();
                    double[] a = Helper.archiveGuideA(seq, sArchive, f, 3);

                    for (int x = 0; x < 4; x++) {
                        pVelocities[i][j][x] = w * pVelocities[i][j][x] +
                                r1 * c1 * (pPersonalBests[i][j][x] - pPositions[i][j][x]) +
                                l * r2 * c2 * (gBest.get(i).getFirst()[x] - pPositions[i][j][x]) +
                                (1 - l) * r3 * c3 * (a[x] - pPositions[i][j][x]);

                        if (vmaxiterlimit > iter) {
                            if (pVelocities[i][j][x] > vmax)
                                pVelocities[i][j][x] = vmax;
                            else if (pVelocities[i][j][x] < -vmax)
                                pVelocities[i][j][x] = -vmax;
                        }

                        pPositions[i][j][x] += pVelocities[i][j][x];
                    }

                    int[][] bitmatrix = Helper.genBitMatrix(pPositions[i][j], seq, colLength);
                    double fitness = f[i].calculate(bitmatrix, seq);
                    if (fitness < pFitnesses[i][j]) {
                        pBitStrings[i][j] = bitmatrixCloner.deepClone(bitmatrix);
                        pPersonalBests[i][j] = positionCloner.deepClone(pPositions[i][j]);
                        pFitnesses[i][j] = fitness;
                    }
                }
            }

            l += 1 / (double) maxIter;

            iter++;
        }
    }

    public static void perform(String[] seq, int[] n, double w, double c1, double c2, double c3, Operator[] ops, int seqNum) throws Exception {
        ArrayList<MGAMPSO> mgampsos = new ArrayList<>();
        FitnessFunction numOfAligned = (bitmatrix, tmpSeq) -> -Helper.numOfAlignedChars(Helper.bitsToStrings(bitmatrix, tmpSeq));
        FitnessFunction insertedIndels = Helper::numOfInsertedIndels;
        FitnessFunction[] f = {numOfAligned, insertedIndels};

        for (int i = 0; i < 30; i++) {
            MGAMPSO ampso = new MGAMPSO(seq, n, w, c1, c2, c3, Double.MAX_VALUE, 0, new double[]{-Double.MAX_VALUE, -Double.MAX_VALUE}, 2500, f, ops);
            ampso.start();
            mgampsos.add(ampso);
        }

        double[][] gBestResultsList = new double[2][30];
        double[] highestResult = {Double.MIN_VALUE, Double.MIN_VALUE};
        double[] lowestResult = {Double.MAX_VALUE, Double.MAX_VALUE};

        int[] infeasibleAmountsList = new int[30];
        int highestInfeasible = Integer.MIN_VALUE;
        int lowestInfeasible = Integer.MAX_VALUE;

        int[][] alignedCharsList = new int[2][30];
        int[] highestAlignedChars = {Integer.MIN_VALUE, Integer.MIN_VALUE};
        int[] lowestAlignedChars = {Integer.MAX_VALUE, Integer.MAX_VALUE};

        int[][] insertedIndelsList = new int[2][30];
        int[] highestInsertedIndels = {Integer.MIN_VALUE, Integer.MIN_VALUE};
        int[] lowestInsertedIndels = {Integer.MAX_VALUE, Integer.MAX_VALUE};

        String[][] bestResultStrings = new String[2][];

        ArrayList<ArrayList<Pair<Integer, Integer>>> pairs = new ArrayList<>();

        for (int i = 0; i < 30; i++) {
            MGAMPSO a = mgampsos.get(i);
            a.join();

            pairs.add(new ArrayList<>());

            Helper.calculateCrowdingDistancesA(seq, a.sArchive, f);
            System.out.println("archive output " + i);
            for (Triplet<double[], int[][], Double> v : a.sArchive) {
                int f1 = Helper.numOfAlignedChars(Helper.bitsToStrings(v.getSecond(), seq));
                int f2 = Helper.numOfInsertedIndels(v.getSecond(), seq);
                System.out.println("(" + f1 + "," + f2 + ") " + Arrays.toString(v.getFirst()) + " " + Arrays.deepToString(v.getSecond()) + " " + v.getThird());
                pairs.get(i).add(new Pair<>(f1, f2));
            }

            double result1 = f[0].calculate(a.gBest.get(0).getSecond(), seq);
            double result2 = f[1].calculate(a.gBest.get(1).getSecond(), seq);
            gBestResultsList[0][i] = result1;
            gBestResultsList[1][i] = result2;

            if (result1 > highestResult[0])
                highestResult[0] = result1;

            if (result2 > highestResult[1])
                highestResult[1] = result2;

            if (result1 < lowestResult[0]) {
                bestResultStrings[0] = Helper.bitsToStrings(a.gBest.get(0).getSecond(), seq);
                lowestResult[0] = result1;
            }

            if (result2 < lowestResult[1]) {
                bestResultStrings[1] = Helper.bitsToStrings(a.gBest.get(1).getSecond(), seq);
                lowestResult[1] = result2;
            }

            int numOfInfeasibleSols = a.numOfInfeasibleSols;
            infeasibleAmountsList[i] = numOfInfeasibleSols;

            if (numOfInfeasibleSols > highestInfeasible)
                highestInfeasible = numOfInfeasibleSols;

            if (numOfInfeasibleSols < lowestInfeasible)
                lowestInfeasible = numOfInfeasibleSols;

            int numOfAlignedChars1 = Helper.numOfAlignedChars(Helper.bitsToStrings(a.gBest.get(0).getSecond(), seq));
            int numOfAlignedChars2 = Helper.numOfAlignedChars(Helper.bitsToStrings(a.gBest.get(1).getSecond(), seq));
            alignedCharsList[0][i] = numOfAlignedChars1;
            alignedCharsList[1][i] = numOfAlignedChars2;

            if (numOfAlignedChars1 > highestAlignedChars[0])
                highestAlignedChars[0] = numOfAlignedChars1;

            if (numOfAlignedChars2 > highestAlignedChars[1])
                highestAlignedChars[1] = numOfAlignedChars2;

            if (numOfAlignedChars1 < lowestAlignedChars[0])
                lowestAlignedChars[0] = numOfAlignedChars1;

            if (numOfAlignedChars2 < lowestAlignedChars[1])
                lowestAlignedChars[1] = numOfAlignedChars2;

            int numOfInsertedIndels1 = Helper.numOfInsertedIndels(a.gBest.get(0).getSecond(), seq);
            int numOfInsertedIndels2 = Helper.numOfInsertedIndels(a.gBest.get(1).getSecond(), seq);
            insertedIndelsList[0][i] = numOfInsertedIndels1;
            insertedIndelsList[1][i] = numOfInsertedIndels2;

            if (numOfInsertedIndels1 > highestInsertedIndels[0])
                highestInsertedIndels[0] = numOfInsertedIndels1;

            if (numOfInsertedIndels2 > highestInsertedIndels[1])
                highestInsertedIndels[1] = numOfInsertedIndels2;

            if (numOfInsertedIndels1 < lowestInsertedIndels[0])
                lowestInsertedIndels[0] = numOfInsertedIndels1;

            if (numOfInsertedIndels2 < lowestInsertedIndels[1])
                lowestInsertedIndels[1] = numOfInsertedIndels2;
        }

        FileOutputStream fos = new FileOutputStream("mgampso-" + seqNum + "-" + n[0] + "-" + w);
        ObjectOutputStream out = new ObjectOutputStream(fos);
        out.writeObject(pairs);

        System.out.println("Final Results:");

        System.out.println("Global Bests:");

        for (int i = 0; i < bestResultStrings.length; i++) {
            System.out.println("[" + i + "]");
            for (String s : bestResultStrings[i]) {
                System.out.println(s);
            }
        }

        System.out.println();

        Helper.printInfo("Fitness", highestResult, lowestResult, gBestResultsList);
        Helper.printInfo("Infeasible Solutions", highestInfeasible, lowestInfeasible, infeasibleAmountsList);
        Helper.printInfo("Alignment", highestAlignedChars, lowestAlignedChars, alignedCharsList);
        Helper.printInfo("Inserted Indels", highestInsertedIndels, lowestInsertedIndels, insertedIndelsList);
    }
}
