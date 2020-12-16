package util;

import java.util.ArrayList;

public class BitMatrix<T extends Number> {
    private ArrayList<ArrayList<Integer>> matrice;

    public BitMatrix(ArrayList<ArrayList<Integer>> matrice) {
        this.matrice = matrice;
    }

    public ArrayList<ArrayList<Integer>> getMatrice() {
        return matrice;
    }
}
