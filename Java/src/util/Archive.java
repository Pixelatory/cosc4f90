package util;

import java.util.ArrayList;

public class Archive<E extends FixedCollection> extends ArrayList<ArrayList<E>> {
    private final int n;

    public Archive(int n, Pair<FitnessFunction, Operator>[] f) {
        this.n = n;
    }

    public void add()
}
