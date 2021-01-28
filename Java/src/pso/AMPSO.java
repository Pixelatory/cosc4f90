package pso;

import base.AM;
import org.kamranzafar.commons.cloner.ObjectCloner;
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
    private double gBestFitness = Double.MIN_VALUE;

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
        else if (term == null || term.length == 0)
            throw new IllegalArgumentException("Termination criteria is empty");
        else if (term.length > 1)
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

        // Initializing the particles of swarm
        for (int i = 0; i < n; i++) {
            int[][] bitstring;
            double[] tmpPos = new double[4];
            double[] tmpVel = new double[4];

            tmpVel[0] = 0;
            tmpVel[1] = 0;
            tmpVel[2] = 0;
            tmpVel[3] = 0;

            do {
                tmpPos[0] = ThreadLocalRandom.current().nextDouble(-0.3, 0.3);
                tmpPos[1] = ThreadLocalRandom.current().nextDouble(0.5, 10);
                tmpPos[2] = ThreadLocalRandom.current().nextDouble(0.5, 10);
                tmpPos[3] = ThreadLocalRandom.current().nextDouble(-0.9, -0.5);

                bitstring = Helper.genBitMatrix(tmpPos, seq, colLength);
            } while (Helper.infeasible(bitstring, seq, ops));

            ObjectCloner<double[]> cloner = new ObjectCloner<>();
            pPositions[i] = tmpPos;
            pPersonalBests[i] = cloner.deepClone(tmpPos);
            pVelocities[i] = tmpVel;
            pBitStrings[i] = bitstring;
            pFitnesses[i] = fitness(bitstring);
        }

        // This is where the iterations begin
        int iter = 0;
        while (iter < maxIter && gBestFitness < term[0]) {
            // Update global best if possible
            for (int i = 0; i < n; i++) {
                if (Helper.infeasible(pBitStrings[i], seq, ops) || (pPositions[i][1] * pPositions[i][2]) == 0)
                    numOfInfeasibleSols++;
                else if (pFitnesses[i] > gBestFitness) {
                    ObjectCloner<int[][]> cloner = new ObjectCloner<>();
                    ObjectCloner<double[]> cloner2 = new ObjectCloner<>();
                    gBestPos = cloner2.deepClone(pPersonalBests[i]);
                    gBestBitString = cloner.deepClone(pBitStrings[i]);
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

                if (!Helper.infeasible(bitstring, seq, ops) && (pPositions[i][1] * pPositions[i][2]) != 0) {
                    double tmpScore = fitness(bitstring);

                    // Update personal best if possible
                    if (tmpScore > pFitnesses[i]) {
                        ObjectCloner<double[]> cloner = new ObjectCloner<>();
                        pPersonalBests[i] = cloner.deepClone(pPositions[i]);
                        pBitStrings[i] = bitstring;
                        pFitnesses[i] = tmpScore;
                    }
                }
            }
            iter++;
        }
    }

    public static void main(String[] args) throws InterruptedException {

        Operator[] ops = {Operator.lt};

        System.out.println("BASIC 1");
        perform(Sequences.basic1, ops);

        System.out.println("BASIC 2");
        perform(Sequences.basic2, ops);

        System.out.println("MED 2");
        perform(Sequences.med2, ops);

        System.out.println("MED 3");
        perform(Sequences.med3, ops);

        System.out.println("SPACES");
        perform(Sequences.spaces, ops);
    }

    public static void perform(String[] seq, Operator[] ops) throws InterruptedException {
        ArrayList<AMPSO> ams = new ArrayList<>();
        int n = 30;
        int maxIter = 5000;
        double w1 = 0.5;
        double w2 = 0.5;

        for (int i = 0; i < 30; i++) {
            AMPSO a = new AMPSO(seq, n, 0.7, 2, 2, Double.MAX_VALUE, 0, new double[]{Double.MAX_VALUE}, maxIter, w1, w2, ops);
            a.start();
            ams.add(a);
        }

        double[] resultList = new double[30];
        double highestResult = Double.MIN_VALUE;
        double lowestResult = Double.MAX_VALUE;

        int[] infeasibleAmountsList = new int[30];
        int highestInfeasible = Integer.MIN_VALUE;
        int lowestInfeasible = Integer.MAX_VALUE;

        int[] alignedCharsList = new int[30];
        int highestAlignedChar = Integer.MIN_VALUE;
        int lowestAlignedChar = Integer.MAX_VALUE;

        int[] insertedIndelsList = new int[30];
        int lowestInsertedIndels = Integer.MAX_VALUE;
        int highestInsertedIndels = Integer.MIN_VALUE;

        String[] bestString = null;

        // Collection information from all of the 30 runs here
        for (int i = 0; i < 30; i++) {
            ams.get(i).join();

            // Overall final fitness
            double result = Helper.aggregatedFunction(ams.get(i).gBestBitString, seq, w1, w2, true, ops);
            resultList[i] = result;

            if (result > highestResult) {
                highestResult = result;
                bestString = Helper.bitsToStrings(ams.get(i).gBestBitString, seq);
            }

            if (result < lowestResult)
                lowestResult = result;

            // Overall number of infeasible solutions
            int numOfInfeasibleSols = ams.get(i).numOfInfeasibleSols;
            infeasibleAmountsList[i] = numOfInfeasibleSols;

            if (numOfInfeasibleSols > highestInfeasible)
                highestInfeasible = numOfInfeasibleSols;

            if (numOfInfeasibleSols < lowestInfeasible)
                lowestInfeasible = numOfInfeasibleSols;

            // Final number of aligned characters
            int numOfAlignedChars = Helper.numOfAlignedChars(Helper.bitsToStrings(ams.get(i).gBestBitString, seq));
            alignedCharsList[i] = numOfAlignedChars;

            if (numOfAlignedChars > highestAlignedChar)
                highestAlignedChar = numOfAlignedChars;

            if (numOfAlignedChars < lowestAlignedChar)
                lowestAlignedChar = numOfAlignedChars;

            // Final number of inserted indels
            int numOfInsertedIndels = Helper.numOfInsertedIndels(ams.get(i).gBestBitString, seq);
            insertedIndelsList[i] = numOfInsertedIndels;

            if (numOfInsertedIndels > highestInsertedIndels)
                highestInsertedIndels = numOfInsertedIndels;

            if (numOfInsertedIndels < lowestInsertedIndels)
                lowestInsertedIndels = numOfInsertedIndels;
        }

        System.out.println("Final Results:");
        System.out.println();

        if (bestString != null) {
            System.out.println("Best result by fitness:");
            for (int i = 0; i < bestString.length; i++) {
                System.out.println(bestString[i]);
            }
        }

        printInfo("Fitness", highestResult, lowestResult, resultList);
        printInfo("Infeasible Solutions", highestInfeasible, lowestInfeasible, infeasibleAmountsList);
        printInfo("Alignment", highestAlignedChar, lowestAlignedChar, alignedCharsList);
        printInfo("Inserted Indels", highestInsertedIndels, lowestInsertedIndels, insertedIndelsList);
    }

    private static void printInfo(String title, double highest, double lowest, double[] list) {
        System.out.println("Best " + title + ": " + highest);
        System.out.println("Worst " + title + ": " + lowest);
        System.out.println("Average " + title + ": " + Helper.average(list));
        System.out.println("St. Dev. " + title + ": " + Helper.stdev(list));
        System.out.println();
    }

    private static void printInfo(String title, int highest, int lowest, int[] list) {
        System.out.println("Highest " + title + ": " + highest);
        System.out.println("Lowest " + title + ": " + lowest);
        System.out.println("Average " + title + ": " + Helper.average(list));
        System.out.println("St. Dev. " + title + ": " + Helper.stdev(list));
        System.out.println();
    }
}
