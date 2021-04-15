package pso;

import base.MGPSO;
import org.kamranzafar.commons.cloner.ObjectCloner;
import util.*;

import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

/*
    Multi-Guided Binary PSO

    Uses linearly increasing lambda value
    from 0-1 over iterations.
 */
public class MGBPSO extends MGPSO {
    private int numOfInfeasibleSols = 0;
    private ArrayList<Pair<int[][], Double>> gBest; // formatted as a pair of bitmatrix and its fitness
    private ArrayList<Pair<int[][], Double>> sArchive;

    public MGBPSO(String[] seq,
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
        startPSO();
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
        int[][][][] pPositions = new int[f.length][highestN][numOfSeqs][colLength];
        int[][][][] pPersonalBests = new int[f.length][highestN][numOfSeqs][colLength];
        double[][] pFitnesses = new double[f.length][highestN];
        double[][][][] pVelocities = new double[f.length][highestN][numOfSeqs][colLength];
        gBest = new ArrayList<>();
        sArchive = new ArrayList<>();
        double l = 0; // lambda coefficient

        ObjectCloner<int[][]> posCloner = new ObjectCloner<>();

        numOfInfeasibleSols = 0;

        // Initialize particles of each sub-swarm
        for (int i = 0; i < f.length; i++) {
            // Within these loops, each particle in each sub-swarm is made
            for (int j = 0; j < n[i]; j++) {

                // This process initializes the particles
                for (int k = 0; k < numOfSeqs; k++) {
                    for (int x = 0; x < colLength; x++) {
                        pPositions[i][j][k][x] = 0;
                        pVelocities[i][j][k][x] = 0;
                    }

                    // This part ensures that all originating particles are feasible
                    for (int y = 0; y < (colLength - seq[k].length()); y++) {
                        while (true) {
                            int randNum = ThreadLocalRandom.current().nextInt(0, colLength);
                            if (pPositions[i][j][k][randNum] != 1) {
                                pPositions[i][j][k][randNum] = 1;
                                break;
                            }
                        }
                    }
                }

                pFitnesses[i][j] = f[i].calculate(pPositions[i][j], seq);
                pPersonalBests[i][j] = posCloner.deepClone(pPositions[i][j]);
            }

            // Just set the global best here, but it'll be changed at the beginning of iterations anyways
            gBest.add(new Pair<>(pPositions[i][0], pFitnesses[i][0]));
        }

        // This is where the iterations begin
        int iter = 0;
        while (iter < maxIter) {
            // If global best is better than termination criteria for respective sub-swarm, then quit
            for (int i = 0; i < f.length; i++) {
                if (gBest.get(i).getSecond() < term[i])
                    return;
            }

            // Update global best and archive if applicable
            for (int i = 0; i < f.length; i++) {
                for (int j = 0; j < n[i]; j++) {
                    // if infeasible, don't try to put into global best
                    if (!Helper.infeasible(pPersonalBests[i][j], seq, ops)) {
                        if (pFitnesses[i][j] < gBest.get(i).getSecond()) {
                            gBest.get(i).setSecond(pFitnesses[i][j]);
                            gBest.get(i).setFirst(posCloner.deepClone(pPersonalBests[i][j]));
                        }
                    }

                    if(!Helper.infeasible(pPositions[i][j], seq, ops)) {
                        Helper.addToArchiveB(seq, sArchive, pPositions[i][j], f, sumN);
                    }
                }
            }

            // Update velocity and position of particles
            for (int i = 0; i < f.length; i++) {
                for (int j = 0; j < n[i]; j++) {
                    double r1 = ThreadLocalRandom.current().nextDouble();
                    double r2 = ThreadLocalRandom.current().nextDouble();
                    double r3 = ThreadLocalRandom.current().nextDouble();
                    int[][] a = Helper.archiveGuideB(seq, sArchive, f, 3);

                    for (int x = 0; x < numOfSeqs; x++) {
                        for (int y = 0; y < colLength; y++) {
                            pVelocities[i][j][x][y] = w * pVelocities[i][j][x][y] +
                                    r1 * c1 * (pPersonalBests[i][j][x][y] - pPositions[i][j][x][y]) +
                                    l * r2 * c2 * (gBest.get(i).getFirst()[x][y] - pPositions[i][j][x][y]) +
                                    (1 - l) * r3 * c3 * (a[x][y] - pPositions[i][j][x][y]);

                            if (vmaxiterlimit > iter) {
                                if (pVelocities[i][j][x][y] > vmax)
                                    pVelocities[i][j][x][y] = vmax;
                                else if (pVelocities[i][j][x][y] < -vmax)
                                    pVelocities[i][j][x][y] = -vmax;
                            }

                            double probability = Helper.Sigmoid(pVelocities[i][j][x][y]);
                            pPositions[i][j][x][y] = ThreadLocalRandom.current().nextDouble() < probability ? 1 : 0;
                        }
                    }

                    // Position is feasible so attempt to update personal best
                    if (!Helper.infeasible(pPositions[i][j], seq, ops)) {
                        double tmpScore = f[i].calculate(pPositions[i][j], seq);
                        if (tmpScore < pFitnesses[i][j]) {
                            pFitnesses[i][j] = tmpScore;
                            pPersonalBests[i][j] = posCloner.deepClone(pPositions[i][j]);
                        }
                    } else
                        numOfInfeasibleSols++;
                }
            }

            // Update lambda parameter (linearly increasing)
            l += 1 / (double) maxIter;

            iter++;
        }
    }

    public static void main(String[] args) throws Exception {
        Operator[] ops = {Operator.lt, Operator.gt};
        String[][] seqss = {Sequences._1bbt_ac, Sequences.yua6_caeel, Sequences.labo_A, Sequences.CSPF_ECOLI, Sequences.SODM_CANAL};
        double[] ws = {0.7548684662307973, 0.664693311042931, 0.6283276486501926, 0.7519652725549882, 0.7499142729146882, 0.7444666103843537, 0.6002649267508061, 0.268990284735547, 0.121495528917823};
        double[] c1s = {1.6403922907668573, 0.12797292579587216, 0.462377069982959, 1.659900433872437, 0.6484731483137582, 1.382159367589628, 1.39119402217411, 1.5805667363521902, 0.49423588598014656};
        double[] c2s = {1.1689889680218306, 1.141633525977906, 0.3400439146349077, 0.17954734318382237, 1.8292436916244121, 0.35852915141793007, 1.7645977300331588, 1.2012673626684758, 0.32362374583425435};
        double[] c3s = {1.0925215379609754, 1.894374658387147, 1.2075173599176188, 0.4740247791781249, 0.47188949238877376, 1.204326868805682, 0.3903042201164242, 1.1353123776616685, 0.5277558146365762};
        int[] n1s = {30, 15, 20, 40};
        int[] n2s = {30, 15, 40, 20};
        for (int i = 0; i < seqss.length; i++) {
            for (int j = 0; j < ws.length; j++) {
                for (int k = 0; k < n1s.length; k++) {
                    System.out.println(ws[j] + " " + c1s[j] + " " + c2s[j] + " " + c3s[j] + " " + n1s[k] + " " + n2s[k]);
                    perform(seqss[i], new int[]{n1s[k], n2s[k]}, ws[j], c1s[j], c2s[j], c3s[j], ops);
                }
            }
        }
    }

    public static void perform(String[] seq, int[] n, double w, double c1, double c2, double c3, Operator[] ops) throws Exception {
        ArrayList<MGBPSO> mgbpsos = new ArrayList<>();
        FitnessFunction numOfAligned = (bitmatrix, tmpSeq) -> -Helper.numOfAlignedChars(Helper.bitsToStrings(bitmatrix, tmpSeq));
        FitnessFunction insertedIndels = Helper::numOfInsertedIndels;
        FitnessFunction[] f = {numOfAligned, insertedIndels};

        for (int i = 0; i < 30; i++) {
            MGBPSO bpso = new MGBPSO(seq, n, w, c1, c2, c3, Double.MAX_VALUE, 0, new double[]{-Double.MAX_VALUE, -Double.MAX_VALUE}, 5000, f, ops);
            bpso.start();
            mgbpsos.add(bpso);
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

        // TODO: Find out what do to with archive

        for (int i = 0; i < 30; i++) {
            MGBPSO b = mgbpsos.get(i);
            b.join();

            double result1 = f[0].calculate(b.gBest.get(0).getFirst(), seq);
            double result2 = f[1].calculate(b.gBest.get(1).getFirst(), seq);
            gBestResultsList[0][i] = result1;
            gBestResultsList[1][i] = result2;

            if (result1 > highestResult[0])
                highestResult[0] = result1;

            if (result2 > highestResult[1])
                highestResult[1] = result2;

            if (result1 < lowestResult[0]) {
                bestResultStrings[0] = Helper.bitsToStrings(b.gBest.get(0).getFirst(), seq);
                lowestResult[0] = result1;
            }

            if (result2 < lowestResult[1]) {
                bestResultStrings[1] = Helper.bitsToStrings(b.gBest.get(1).getFirst(), seq);
                lowestResult[1] = result2;
            }

            int numOfInfeasibleSols = b.numOfInfeasibleSols;
            infeasibleAmountsList[i] = numOfInfeasibleSols;

            if (numOfInfeasibleSols > highestInfeasible)
                highestInfeasible = numOfInfeasibleSols;

            if (numOfInfeasibleSols < lowestInfeasible)
                lowestInfeasible = numOfInfeasibleSols;

            int numOfAlignedChars1 = Helper.numOfAlignedChars(Helper.bitsToStrings(b.gBest.get(0).getFirst(), seq));
            int numOfAlignedChars2 = Helper.numOfAlignedChars(Helper.bitsToStrings(b.gBest.get(1).getFirst(), seq));
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

            int numOfInsertedIndels1 = Helper.numOfInsertedIndels(b.gBest.get(0).getFirst(), seq);
            int numOfInsertedIndels2 = Helper.numOfInsertedIndels(b.gBest.get(1).getFirst(), seq);
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
