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

        numOfInfeasibleSols = 0;

        // Initialize particles of each sub-swarm
        for (int i = 0; i < f.length; i++) {
            // These are the containers for each sub-swarm
            int[][][] newPositions = new int[n[i]][numOfSeqs][colLength];
            double[][][] newVelocities = new double[n[i]][numOfSeqs][colLength];

            // Within these loops, each particle in each sub-swarm is made
            for (int j = 0; j < n[i]; j++) {
                int[][] tmpPos = new int[numOfSeqs][colLength];
                double[][] tmpVel = new double[numOfSeqs][colLength];

                // This process initializes the particles
                for (int k = 0; k < numOfSeqs; k++) {
                    for (int x = 0; x < colLength; x++) {
                        tmpPos[k][x] = 0;
                        tmpVel[k][x] = 0;
                    }

                    // This part ensures that all originating particles are feasible
                    for (int y = 0; y < (colLength - seq[k].length()); y++) {
                        while (true) {
                            int randNum = ThreadLocalRandom.current().nextInt(0, colLength);
                            if (tmpPos[k][randNum] != 1) {
                                tmpPos[k][randNum] = 1;
                                break;
                            }
                        }
                    }
                }

                newPositions[j] = tmpPos;
                newVelocities[j] = tmpVel;
                pFitnesses[i][j] = f[i].calculate(tmpPos, seq);
            }

            ObjectCloner<int[][][]> cloner = new ObjectCloner<>();
            pPositions[i] = newPositions;
            pVelocities[i] = newVelocities;
            pPersonalBests[i] = cloner.deepClone(newPositions);

            for (int[][] position : newPositions) {
                if (gBest.size() < (i + 1)) {
                    gBest.add(new Pair<>(position, f[i].calculate(position, seq)));
                } else {
                    double tmpScore = f[i].calculate(position, seq);
                    if (tmpScore > gBest.get(i).getSecond()) {
                        ObjectCloner<int[][]> cloner2 = new ObjectCloner<>();
                        gBest.set(i, new Pair<>(cloner2.deepClone(position), tmpScore));
                    }
                }
            }
        }

        // This is where the iterations begin
        int iter = 0;
        while (iter < maxIter) {
            // If global best is better than termination criteria for respective sub-swarm, then quit
            for (int i = 0; i < f.length; i++) {
                if (gBest.get(i).getSecond() < term[i])
                    return;
            }

            // Update global best if applicable
            for (int i = 0; i < f.length; i++) {
                for (int j = 0; j < n[i]; j++) {
                    Helper.addToArchive(seq, sArchive, pPositions[i][j], f, sumN);

                    // if infeasible, don't try to put into global best
                    if (Helper.infeasible(pPositions[i][j], seq, ops)) {
                        numOfInfeasibleSols++;
                        continue;
                    }

                    double tmpScore = f[i].calculate(pPositions[i][j], seq);

                    if (tmpScore > gBest.get(i).getSecond()) {
                        ObjectCloner<int[][]> cloner = new ObjectCloner<>();
                        gBest.get(i).setSecond(tmpScore);
                        gBest.get(i).setFirst(cloner.deepClone(pPersonalBests[i][j]));
                    }
                }
            }

            // Update velocity and position of particles
            for (int i = 0; i < f.length; i++) {
                for (int j = 0; j < n[i]; j++) {
                    double r1 = ThreadLocalRandom.current().nextDouble();
                    double r2 = ThreadLocalRandom.current().nextDouble();
                    double r3 = ThreadLocalRandom.current().nextDouble();
                    int[][] a = Helper.archiveGuide(seq, sArchive, f, 3);

                    for (int x = 0; x < numOfSeqs; x++) {
                        for (int y = 0; y < colLength; y++) {
                            pVelocities[i][j][x][y] = w * pVelocities[i][j][x][y] +
                                    r1 * c1 * (pPersonalBests[i][j][x][y] - pPositions[i][j][x][y]) +
                                    l * r2 * c2 * (gBest.get(i).getFirst()[x][y] - pPositions[i][j][x][y]) +
                                    (1 - l) * r3 * c3 * (a[x][y] - pPositions[i][j][x][y]);

                            if (vmaxiterlimit < iter) {
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
                        if (tmpScore > pFitnesses[i][j]) {
                            ObjectCloner<int[][]> cloner = new ObjectCloner<>();
                            pFitnesses[i][j] = tmpScore;
                            pPersonalBests[i][j] = cloner.deepClone(pPositions[i][j]);
                        }
                    }
                }
            }

            // Update lambda parameter (linearly increasing)
            l += 1 / (double) maxIter;

            iter++;
        }
    }

    public static void main(String[] args) throws Exception {
        FitnessFunction numOfAligned = (bitmatrix, seq) -> -Helper.numOfAlignedChars(Helper.bitsToStrings(bitmatrix, seq));
        FitnessFunction insertedIndels = Helper::numOfInsertedIndels;
        MGBPSO mg = new MGBPSO(Sequences.basic1, new int[]{20, 30}, 0.75, 1.0, 1.6, 1.05, Double.MAX_VALUE, Integer.MAX_VALUE, new double[]{-Double.MAX_VALUE, -Double.MAX_VALUE}, 5000, new FitnessFunction[]{numOfAligned, insertedIndels}, new Operator[]{Operator.lt});
        mg.startPSO();

        System.out.println("Global Bests:");
        for (Pair<int[][], Double> m : mg.gBest) {
            System.out.print(Math.abs(numOfAligned.calculate(m.getFirst(), Sequences.basic1)));
            System.out.print(" " + insertedIndels.calculate(m.getFirst(), Sequences.basic1));
            System.out.println();
            for (String s : Helper.bitsToStrings(m.getFirst(), Sequences.basic1)) {
                System.out.println(s);
            }
        }

        System.out.println("ARCHIVE");
        for (Pair<int[][], Double> s : mg.sArchive) {
            System.out.print(Math.abs(numOfAligned.calculate(s.getFirst(), Sequences.basic1)));
            System.out.print(" " + insertedIndels.calculate(s.getFirst(), Sequences.basic1));
            System.out.println();

            for (String ss : Helper.bitsToStrings(s.getFirst(), Sequences.basic1))
                System.out.println(ss);
        }
    }
}
