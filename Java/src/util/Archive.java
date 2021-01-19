/*package util;

import java.util.ArrayList;

public class Archive<E extends FixedCollection> {
    private final int n;
    private ArrayList<E> archive;
    private final Pair<FitnessFunction, Operator>[] f;

    public Archive(int n, Pair<FitnessFunction, Operator>[] f) {
        this.n = n;
        this.archive = new ArrayList<>();
        this.f = f;
    }

    /**
     * Checks if one bitmatrix dominates another.
     * @param seq
     * @param bm1
     * @param bm2
     * @return
     */
/*
    private static boolean dominates(ArrayList<String> seq,
                                     ArrayList<ArrayList<Integer>> bm1,
                                     ArrayList<ArrayList<Integer>> bm2) {
        ArrayList<String> strings1 = Helper.bitsToStrings(bm1, seq);
        ArrayList<String> strings2 = Helper.bitsToStrings(bm2, seq);

        // First check the number of aligned characters
        int res1 = Helper.numOfAlignedChars(strings1);
        int res2 = Helper.numOfAlignedChars(strings2);

        if(res1 <= res2)
            return false;

        // The code won't reach here if bm1 already doesn't have more aligned characters than bm2
        // Now check the number of inserted indels
        res1 = Helper.numOfInsertedIndels(bm1, seq);
        res2 = Helper.numOfInsertedIndels(bm2, seq);

        return res1 > res2;
    }

    public void add(int fIndex, E e) {
        for(E entry : archive) {
            if()
        }
        archive.get(fIndex).add(e);
    }
}
*/