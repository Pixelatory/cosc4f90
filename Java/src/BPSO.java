import util.Helper;
import util.Operator;
import util.PSO;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.concurrent.ThreadLocalRandom;

public class BPSO extends PSO {
    public BPSO(ArrayList<String> seq,
                int n,
                double w,
                double c1,
                double c2,
                double vmax,
                int vmaxiterlimit,
                double term,
                int maxIter,
                ArrayList<Operator> ops) {
        super(seq, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, ops);
    }

    public void start() {
        // Begin checks for trivial errors
        if (n < 1)
            throw new IllegalArgumentException("Swarm size cannot be < 1");
        else if (seq.size() < 2)
            throw new IllegalArgumentException("Number of sequences < 2");
        else if (maxIter < 1)
            throw new IllegalArgumentException("Maximum iterations cannot be < 1");

        // Initialize main data containers
        ArrayList<ArrayList<ArrayList<Integer>>> gBestPos = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Integer>>> pPositions = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Integer>>> pPersonalBests = new ArrayList<>();
        ArrayList<ArrayList<ArrayList<Double>>> pVelocities = new ArrayList<>();
        ArrayList<Double> pFitnesses = new ArrayList<>(); // so we only calculate fitness once
        int numOfInfeasibleSols = 0;
        int gBestFitness = 0;

        // More vars just to make things more readable
        int colLength = Helper.getColLength(seq);
        int numOfSeq = seq.size();

        // Initializing the particles of the swarm
        for (int i = 0; i < n; i++) {
            ArrayList<ArrayList<Integer>> position = new ArrayList<>();
            ArrayList<ArrayList<Double>> velocity = new ArrayList<>();

            for (int j = 0; j < numOfSeq; j++) {
                velocity.add(new ArrayList<>(colLength));
                position.add(new ArrayList<>(colLength));
                for (int k = 0; k < colLength; k++) {
                    velocity.get(j).add(0.0);
                    position.get(j).add(0);
                }
                for (int k = 0; k < colLength - seq.get(j).length(); k++) {
                    while(true) {
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
        }
    }

    public static void main(String[] args) {
        ArrayList<String> seqs = new ArrayList<>();
        seqs.add("22222");
        seqs.add("2515125");
        BPSO b = new BPSO(seqs, 30, 1.4, 1.4, 1.4, 11, 100, 100, 100, null);
        b.start();
    }
}
