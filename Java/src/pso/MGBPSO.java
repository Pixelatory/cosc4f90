package pso;

import base.MGPSO;
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
                  int n,
                  double w,
                  double c1,
                  double c2,
                  double c3,
                  double vmax,
                  int vmaxiterlimit,
                  double[] term,
                  int maxIter,
                  FitnessFunction[] f,
                  Operator[] ops) {
        super(seq, n, w, c1, c2, c3, vmax, vmaxiterlimit, term, maxIter, f, ops);
    }

    public void startPSO() {
        // Begin checks for trivial errors
        if (n < 1)
            throw new IllegalArgumentException("Swarm size cannot be < 1");
        else if (seq.length < 2)
            throw new IllegalArgumentException("Number of sequences < 2");
        else if (maxIter < 1)
            throw new IllegalArgumentException("Maximum iterations cannot be < 1");

        int colLength = Helper.getColLength(seq);
        int numOfSeqs = seq.length;

        // Initialize main data containers
        int[][][][] pPositions = new int[f.size()][n][numOfSeqs][colLength];
        int[][][][] pPersonalBests = new int[f.size()][n][numOfSeqs][colLength];
        double[][] pFitnesses = new double[f.size()][n];
        double[][][][] pVelocities = new double[f.size()][n][numOfSeqs][colLength];
        gBest = new ArrayList<>();
        sArchive = new ArrayList<>();
        double l = 0; // lambda coefficient

        numOfInfeasibleSols = 0;

        // Initialize the global best (one per swarm, or per entry in f)
        for (int i = 0; i < f.size(); i++) {
            int[][] tmpPos = new int[numOfSeqs][colLength];

            // gBest Position starts at all 0s
            for (int j = 0; j < numOfSeqs; j++) {
                for (int k = 0; k < colLength; k++) {
                    tmpPos[i][j] = 0;
                }
            }

            gBest.add(new Pair<>(tmpPos, f.get(i).getFirst().calculate(tmpPos, seq)));
        }

        // Initialize particles of each sub-swarm
        for (int i = 0; i < f.size(); i++) {
            // These are the containers for each sub-swarm
            int[][][] newPositions = new int[n][numOfSeqs][colLength];
            double[][][] newVelocities = new double[n][numOfSeqs][colLength];

            // Within these loops, each particle in each sub-swarm is made
            for (int j = 0; j < n; j++) {
                int[][] tmpPos = new int[numOfSeqs][colLength];
                double[][] tmpVel = new double[numOfSeqs][colLength];

                // This process initializes the particles
                for (int k = 0; k < numOfSeqs; k++) {
                    for (int x = 0; x < colLength; x++) {
                        tmpPos[k][x] = 0;
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
                pFitnesses[i][j] = f.get(i).getFirst().calculate(tmpPos, seq);
            }

            pPositions[i] = newPositions;
            pVelocities[i] = newVelocities;
            pPersonalBests[i] = ArrayCloner.deepcopy(newPositions);
        }

        // This is where the iterations begin
        int iter = 0;
        while (iter < maxIter) {
            // If global best is better than termination criteria for respective sub-swarm, then quit
            for (int i = 0; i < f.size(); i++) {
                if (gBest.get(i).getSecond() < term[i])
                    return;
            }

            // Update global best if applicable
            for (int i = 0; i < f.size(); i++) {
                for (int j = 0; j < n; j++) {
                    // if infeasible, don't try to put into global best or archive
                    if (Helper.infeasible(pPositions[i][j], seq, ops)) {
                        numOfInfeasibleSols++;
                        continue;
                    }

                    sArchive = Helper.addToArchive(seq, sArchive, pPositions[i][j], f);

                    if (!Helper.infeasible(pPersonalBests[i][j], seq, ops)
                            && pFitnesses[i][j] > gBest.get(i).getSecond()) {
                        gBest.get(i).setSecond(pFitnesses[i][j]);
                        gBest.get(i).setFirst(ArrayCloner.deepcopy(pPersonalBests[i][j]));
                    }
                }
            }

            for (int i = 0; i < f.size(); i++) {
                for (int j = 0; j < n; j++) {
                    double r1 = ThreadLocalRandom.current().nextDouble();
                    double r2 = ThreadLocalRandom.current().nextDouble();
                    double r3 = ThreadLocalRandom.current().nextDouble();

                }
            }

            // Update lambda parameter (linearly increasing)
            l += 1 / (double) maxIter;

            iter++;
        }
    }
}
