import util.Helper;
import util.Operator;
import base.PSO;

import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

public class BPSO extends PSO {
    private final double w1;
    private final double w2;
    private int numOfInfeasibleSols = 0;
    private ArrayList<ArrayList<Integer>> gBestPos;

    private double fitness(ArrayList<ArrayList<Integer>> bitmatrix) {
        return Helper.aggregatedFunction(bitmatrix, seq, w1, w2, true, ops);
    }

    public BPSO(ArrayList<String> seq,
                int n,
                double w,
                double c1,
                double c2,
                double vmax,
                int vmaxiterlimit,
                double term,
                int maxIter,
                double w1,
                double w2,
                ArrayList<Operator> ops) {
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
        else if (seq.size() < 2)
            throw new IllegalArgumentException("Number of sequences < 2");
        else if (maxIter < 1)
            throw new IllegalArgumentException("Maximum iterations cannot be < 1");

        // Initialize main data containers
        gBestPos = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Integer>>> pPositions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Integer>>> pPersonalBests = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> pVelocities = new ArrayList<>();
        ArrayList<Double> pFitnesses = new ArrayList<>(); // so we only calculate fitness once
        numOfInfeasibleSols = 0;
        double gBestFitness;

        // More vars just to make things more readable
        int colLength = Helper.getColLength(seq);
        int numOfSeq = seq.size();

        // Initializing the global best particle to all 0 integer bitmatrix
        for (int i = 0; i < numOfSeq; i++) {
            ArrayList<Integer> tmp = new ArrayList<>();
            for (int j = 0; j < colLength; j++) {
                tmp.add(0);
            }
            gBestPos.add(tmp);
        }

        gBestFitness = fitness(gBestPos); // Initializing the global best fitness

        // Initializing the particles of the swarm (they will always initially be feasible)
        for (int i = 0; i < n; i++) {
            ArrayList<ArrayList<Integer>> position = new ArrayList<>();
            ArrayList<ArrayList<Double>> velocity = new ArrayList<>();

            for (int j = 0; j < numOfSeq; j++) {
                velocity.add(new ArrayList<>());
                position.add(new ArrayList<>());
                for (int k = 0; k < colLength; k++) {
                    velocity.get(j).add(0.0);
                    position.get(j).add(0);
                }
                for (int k = 0; k < colLength - seq.get(j).length(); k++) {
                    while (true) {
                        int randNum = ThreadLocalRandom.current().nextInt(0, colLength);
                        if (position.get(j).get(randNum) != 1) {
                            position.get(j).set(randNum, 1);
                            break;
                        }
                    }
                }
            }

            pPositions.add(position);
            pPersonalBests.add(position);
            pVelocities.add(velocity);
            double tmpParticleFitness = fitness(position);
            pFitnesses.add(tmpParticleFitness);

            if (tmpParticleFitness > gBestFitness) {
                gBestPos = Helper.copyArray(position);
                gBestFitness = tmpParticleFitness;
            }
        }

        // This is where the iterations begin
        int iter = 0; // iteration counter
        while (iter < maxIter && fitness(gBestPos) < term) {

            // Checking of fitness is better than global best for updating
            for (int i = 0; i < n; i++) {
                if (pFitnesses.get(i) > gBestFitness) {
                    gBestPos = Helper.copyArray(pPersonalBests.get(i));
                    gBestFitness = pFitnesses.get(i);
                }
            }

            // Update each particle's velocity, position, and personal best
            for (int i = 0; i < n; i++) {
                // r1 and r2 are ~ U(0,1)
                double r1 = ThreadLocalRandom.current().nextDouble(0, 1);
                double r2 = ThreadLocalRandom.current().nextDouble(0, 1);

                for (int j = 0; j < numOfSeq; j++) {
                    for (int k = 0; k < colLength; k++) {
                        pVelocities.get(i).get(j).set(k, w * pVelocities.get(i).get(j).get(k)
                                + r1 * c1 * (pPersonalBests.get(i).get(j).get(k) - pPositions.get(i).get(j).get(k))
                                + r2 * c2 * (gBestPos.get(j).get(k) - pPositions.get(i).get(j).get(k)));

                        if (vmaxiterlimit < iter) {
                            if (pVelocities.get(i).get(j).get(k) > vmax)
                                pVelocities.get(i).get(j).set(k, vmax);
                            else if (pVelocities.get(i).get(j).get(k) < -vmax)
                                pVelocities.get(i).get(j).set(k, -vmax);
                        }

                        double probability = Helper.Sigmoid(pVelocities.get(i).get(j).get(k));
                        pPositions.get(i).get(j).set(k, ThreadLocalRandom.current().nextDouble(0, 1) < probability ? 1 : 0);
                    }
                }

                double tmpFitness = fitness(pPositions.get(i));

                // solution is infeasible, so increment count
                if (tmpFitness == Double.MIN_VALUE) {
                    numOfInfeasibleSols++;
                    continue;
                }

                // Checking if the fitness is better than personal best for updating
                if (tmpFitness > pFitnesses.get(i)) {
                    pFitnesses.set(i, tmpFitness);
                    pPersonalBests.set(i, Helper.copyArray(pPositions.get(i)));
                }
            }

            iter++;
        }
    }

    public static void main(String[] args) throws InterruptedException {
        ArrayList<BPSO> bs = new ArrayList<>();

        Operator lt = (a, b) -> a < b;
        ArrayList<Operator> ops = new ArrayList<>();
        ops.add(lt);

        ArrayList<String> seqs = new ArrayList<>();
        seqs.add("CBCADCAACE");
        seqs.add("EACABDCADB");
        seqs.add("DABAECBDCD");
        seqs.add("DBEACEACCD");
        seqs.add("DDABDEEEDE");
        seqs.add("EEAECCAAEB");
        seqs.add("EABEBCBCCB");
        seqs.add("BAADDACDBB");

        int n = 30;
        int maxIter = 5000;
        double w1 = 0.5;
        double w2 = 0.5;

        for (int i = 0; i < 30; i++) {
            BPSO b = new BPSO(seqs, n, 0.99, 2, 2, 11, 0, Double.MAX_VALUE, maxIter, w1, w2, ops);
            b.start();
            bs.add(b);
        }

        for (int i = 0; i < 30; i++) {
            bs.get(i).join();
            System.out.println("Test Run:");
            System.out.println((((double) bs.get(i).numOfInfeasibleSols / (n * maxIter)) * 100) + "% infeasible");
            System.out.println(Helper.aggregatedFunction(bs.get(i).gBestPos, seqs, w1, w2, true, ops));
        }


    }
}