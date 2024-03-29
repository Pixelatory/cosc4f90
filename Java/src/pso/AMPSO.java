package pso;

import base.AM;
import org.kamranzafar.commons.cloner.ObjectCloner;
import util.Helper;
import util.Operator;
import util.TFAParser;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;

public class AMPSO extends AM {
    private final double w1;
    private final double w2;
    private int numOfInfeasibleSols = 0;
    private double[] gBestPos;
    private int[][] gBestBitString;
    private double gBestFitness = Double.MIN_VALUE;
    private double[] popDiversities;

    private double fitness(int[][] bitmatrix) {
        //return Helper.aggregatedFunction(bitmatrix, seq, w1, w2, true, ops);
        return Helper.Similarity(bitmatrix, seq, true, ops);
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
        this.popDiversities = new double[maxIter];
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

        // Vars just to make things more readable
        int colLength = Helper.getColLength(seq);
        int numOfSeqs = seq.length;

        // Initialize data containers
        double[][] pPositions = new double[n][4];
        double[][] pPersonalBests = new double[n][4];
        double[][] pVelocities = new double[n][4];
        int[][][] pBitStrings = new int[n][numOfSeqs][colLength];
        double[] pFitnesses = new double[n]; // so we only calculate fitness once

        // Initialize cloners
        ObjectCloner<double[]> positionCloner = new ObjectCloner<>();
        ObjectCloner<int[][]> bitmatrixCloner = new ObjectCloner<>();

        numOfInfeasibleSols = 0;

        // Initializing the particles of swarm
        for (int i = 0; i < n; i++) {
            int[][] bitstring;
            pVelocities[i][0] = 0;
            pVelocities[i][1] = 0;
            pVelocities[i][2] = 0;
            pVelocities[i][3] = 0;

            do {
                pPositions[i][0] = ThreadLocalRandom.current().nextDouble(-1, 1);
                pPositions[i][1] = ThreadLocalRandom.current().nextDouble(-1, 1);
                pPositions[i][2] = ThreadLocalRandom.current().nextDouble(-1, 1);
                pPositions[i][3] = ThreadLocalRandom.current().nextDouble(-1, 1);

                bitstring = Helper.genBitMatrix(pPositions[i], seq, colLength);
            } while (Helper.infeasible(bitstring, seq, ops)
                    || pPositions[i][1] * pPositions[i][2] == 0);

            pPersonalBests[i] = positionCloner.deepClone(pPositions[i]);
            pBitStrings[i] = bitstring;
            pFitnesses[i] = fitness(bitstring);
        }

        // This is where the iterations begin
        int iter = 0;
        while (iter < maxIter && gBestFitness < term[0]) {

            // Updating global best if possible
            /*for (int i = 0; i < n; i++) {
                if (pFitnesses[i] > gBestFitness) {
                    gBestPos = positionCloner.deepClone(pPersonalBests[i]);
                    gBestBitString = bitmatrixCloner.deepClone(pBitStrings[i]);
                    gBestFitness = pFitnesses[i];
                }
            }*/

            // Calculate diversity of swarm positions
            double pPositionDiversity = 0;
            for (int i = 0; i < 4; i++) {
                // Dimension mean
                double sum = 0;
                for (int j = 0; j < n; j++) {
                    sum += pPositions[j][i];
                }

                double dimensionDivMean = sum / n;

                // Population diversity on this dimension
                sum = 0;
                for (int j = 0; j < n; j++) {
                    sum += Math.abs(pPositions[j][i] - dimensionDivMean);
                }
                pPositionDiversity += sum / n;
            }
            // Total population diversity
            popDiversities[iter] = pPositionDiversity / 4;

            // Update each particle's velocity, position, and personal best
            for (int i = 0; i < n; i++) {
                double r1 = ThreadLocalRandom.current().nextDouble();
                double r2 = ThreadLocalRandom.current().nextDouble();

                // Get the indice positions for neighbourhood
                int[] idxs = {i - 1, i, i + 1};

                // The indices wrap around if outside of array boundary.
                if (i - 1 < 0)
                    idxs[0] = n - 1;
                if (i + 1 >= n)
                    idxs[2] = 0;

                double bestFitness = pFitnesses[idxs[0]];
                double[] lBestPosition = pPersonalBests[idxs[0]];
                if (pFitnesses[idxs[1]] > bestFitness) {
                    bestFitness = pFitnesses[idxs[1]];
                    lBestPosition = pPersonalBests[idxs[1]];
                }

                if (pFitnesses[idxs[2]] > bestFitness) {
                    lBestPosition = pPersonalBests[idxs[2]];
                }

                // Velocity update
                for (int j = 0; j < 4; j++) {
                    pVelocities[i][j] = w * pVelocities[i][j] +
                            r1 * c1 * (pPersonalBests[i][j] - pPositions[i][j]) +
                            r2 * c2 * (lBestPosition[j] - pPositions[i][j]);

                    if (vmaxiterlimit > iter) {
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

                if (Helper.infeasible(bitstring, seq, ops)
                        || pPositions[i][1] * pPositions[i][2] == 0)
                    numOfInfeasibleSols++;
                else if (tmpScore > pFitnesses[i]) {
                    pPersonalBests[i] = positionCloner.deepClone(pPositions[i]);
                    pBitStrings[i] = bitstring;
                    pFitnesses[i] = tmpScore;
                }
            }
            iter++;
        }

        // Iterations complete, so set the global best result here
        for (int i = 0; i < pFitnesses.length; i++) {
            if (pFitnesses[i] > gBestFitness) {
                gBestFitness = pFitnesses[i];
                gBestPos = pPersonalBests[i];
                gBestBitString = pBitStrings[i];
            }
        }
    }

    public static void main(String[] args) throws InterruptedException, IOException {
        System.setOut(new PrintStream("ampso lBest Similarity as fitness with diversity measurement.txt"));
        Operator[] ops = {Operator.lt};
        /*String[][] seqss = {TFAParser.seq1,
                TFAParser.seq2,
                TFAParser.seq3,
                TFAParser.seq4,
                TFAParser.seq5,
                TFAParser.seq6,
                TFAParser.seq7,
                TFAParser.seq8,
                TFAParser.seq9,
                TFAParser.seq10,
                TFAParser.seq11,
                TFAParser.seq12};*/

        String[][] seqss = {TFAParser.OneAMK};
        double[] ws = {0.7098150314023034, 0.7861878289341226, 0.7226988763584911, 0.7566563010222623, 0.3260770959583924};
        double[] c1s = {1.6788775458244407, 1.2489605895810598, 0.7877525185622731, 0.9360167168210383, 1.397180934100446};
        double[] c2s = {0.899446463381824, 1.211314077689869, 1.3569405150479763, 0.8238319031198684, 1.3464223101885027};
        double[] w1s = {1, 0, 0.5};
        double[] w2s = {0, 1, 0.5};

        for (int i = 0; i < seqss.length; i++) {
            for (int j = 0; j < ws.length; j++) {
                //for (int k = 0; k < w1s.length; k++) {
                System.out.println(i + " " + 50 + " " + ws[j] + " " + c1s[j] + " " + c2s[j] + " " + 1 + " " + 0);
                perform(seqss[i], ops, 50, ws[j], c1s[j], c2s[j], 0, 0);
                //}
            }
        }
        //Process child = Runtime.getRuntime().exec("shutdown -s");
    }

    public static void perform(String[] seq, Operator[] ops, int n, double w, double c1, double c2, double w1, double w2) throws InterruptedException {
        ArrayList<AMPSO> ams = new ArrayList<>();
        int maxIter = 2500;

        for (int i = 0; i < 30; i++) {
            AMPSO a = new AMPSO(seq, n, w, c1, c2, Double.MAX_VALUE, 0, new double[]{Double.MAX_VALUE}, maxIter, w1, w2, ops);
            a.start();
            ams.add(a);
            ams.get(i).join();
        }

        double[] avgDiversities = new double[maxIter];
        for (int i = 0; i < maxIter; i++) {
            double sum = 0;
            for (int j = 0; j < 30; j++) {
                sum += ams.get(j).popDiversities[i];
            }
            avgDiversities[i] = sum / 30;
        }

        System.out.println(Arrays.toString(avgDiversities));

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

            // Overall final fitness
            //double result = Helper.aggregatedFunction(ams.get(i).gBestBitString, seq, w1, w2, true, ops);
            double result = Helper.Similarity(ams.get(i).gBestBitString, seq, true, ops);
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

        Helper.printInfo("Fitness", highestResult, lowestResult, resultList);
        Helper.printInfo("Infeasible Solutions", highestInfeasible, lowestInfeasible, infeasibleAmountsList);
        Helper.printInfo("Alignment", highestAlignedChar, lowestAlignedChar, alignedCharsList);
        Helper.printInfo("Inserted Indels", highestInsertedIndels, lowestInsertedIndels, insertedIndelsList);
    }
}
