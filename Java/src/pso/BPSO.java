package pso;

import org.kamranzafar.commons.cloner.ObjectCloner;
import util.*;
import base.PSO;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.IntStream;

public class BPSO extends PSO {
    private final double w1;
    private final double w2;
    private int[][] gBestPos;
    private double gBestFitness;
    private int infeasibleSols; // Holds count of infeasible solutions from each iteration
    private ArrayList<Pair<Integer, Double>> gBests; // Holds global best fitness and iteration value of its superiority

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
        infeasibleSols = 0;

        gBests = new ArrayList<>();

        ObjectCloner<int[][]> posCloner = new ObjectCloner<>();

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

            pPositions[i] = position;
            pPersonalBests[i] = posCloner.deepClone(position);
            pVelocities[i] = velocity;
            pFitnesses[i] = fitness(position);
        }

        // This is where the iterations begin
        int iter = 0; // iteration counter
        while (iter < maxIter && gBestFitness < term[0]) {

            // Checking if fitness is better than global best for updating
            for (int i = 0; i < n; i++) {
                if (pFitnesses[i] > gBestFitness) {
                    gBestPos = posCloner.deepClone(pPersonalBests[i]);
                    gBestFitness = pFitnesses[i];
                    gBests.add(new Pair<>(iter, gBestFitness));
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

                        System.out.println(pVelocities[i][j][k]);

                        if (vmaxiterlimit > iter) {
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
                if (Helper.infeasible(pPositions[i], seq, ops))
                    infeasibleSols++;
                else if (tmpFitness > pFitnesses[i]) {
                    // current fitness is better than personal best's so update it
                    pFitnesses[i] = tmpFitness;
                    pPersonalBests[i] = posCloner.deepClone(pPositions[i]);
                }
            }

            iter++;
        }
    }

    public static void main(String[] args) throws InterruptedException, IOException {
        //System.setOut(new PrintStream("output-file-bpso-ltgt-preliminary.txt"));
        Operator[] ops = {Operator.lt, Operator.gt};
        String[][] seqss = {TFAParser.seq1,
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
                TFAParser.seq12};
        double[] ws = {0.7098150314023034, 0.7861878289341226, 0.7226988763584911, 0.7566563010222623, 0.3260770959583924};
        double[] c1s = {1.6788775458244407, 1.2489605895810598, 0.7877525185622731, 0.9360167168210383, 1.397180934100446};
        double[] c2s = {0.899446463381824, 1.211314077689869, 1.3569405150479763, 0.8238319031198684, 1.3464223101885027};
        double[] w1s = {1, 0, 0.5};
        double[] w2s = {0, 1, 0.5};
        int[] ns = {30, 40, 50};

        for (int i = 0; i < seqss.length; i++) {
            for (int t = 0; t < ns.length; t++) {
                for (int j = 0; j < ws.length; j++) {
                    for (int k = 0; k < w1s.length; k++) {
                        System.out.println(i + " " + ns[t] + " " + ws[j] + " " + c1s[j] + " " + c2s[j] + " " + w1s[k] + " " + w2s[k]);
                        perform(seqss[i], ops, ns[t], ws[j], c1s[j], c2s[j], w1s[k], w2s[k]);
                    }
                }
            }
        }

        Process child = Runtime.getRuntime().exec("shutdown -s");
    }

    public static void perform(String[] seq, Operator[] ops, int n, double w, double c1, double c2, double w1, double w2) throws InterruptedException {
        ArrayList<BPSO> bs = new ArrayList<>();

        int maxIter = 2500;

        for (int i = 0; i < 30; i++) {
            BPSO b = new BPSO(seq, n, w, c1, c2, Double.MAX_VALUE, 0, new double[]{Double.MAX_VALUE}, maxIter, w1, w2, ops);
            b.start();
            bs.add(b);
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
            bs.get(i).join();

            // Overall final fitness
            double result = Helper.aggregatedFunction(bs.get(i).gBestPos, seq, w1, w2, true, ops);
            resultList[i] = result;

            if (result > highestResult) {
                highestResult = result;
                bestString = Helper.bitsToStrings(bs.get(i).gBestPos, seq);
            }

            if (result < lowestResult)
                lowestResult = result;

            // Overall number of infeasible solutions
            int numOfInfeasibleSols = IntStream.of(bs.get(i).infeasibleSols).sum();
            infeasibleAmountsList[i] = numOfInfeasibleSols;

            if (numOfInfeasibleSols > highestInfeasible)
                highestInfeasible = numOfInfeasibleSols;

            if (numOfInfeasibleSols < lowestInfeasible)
                lowestInfeasible = numOfInfeasibleSols;

            // Final number of aligned characters
            int numOfAlignedChars = Helper.numOfAlignedChars(Helper.bitsToStrings(bs.get(i).gBestPos, seq));
            alignedCharsList[i] = numOfAlignedChars;

            if (numOfAlignedChars > highestAlignedChar)
                highestAlignedChar = numOfAlignedChars;

            if (numOfAlignedChars < lowestAlignedChar)
                lowestAlignedChar = numOfAlignedChars;

            // Final number of inserted indels
            int numOfInsertedIndels = Helper.numOfInsertedIndels(bs.get(i).gBestPos, seq);
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