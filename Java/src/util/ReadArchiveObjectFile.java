package util;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;

public class ReadArchiveObjectFile {
    public static void main(String[] args) throws IOException, ClassNotFoundException {
        ObjectInputStream inputStream = new ObjectInputStream(new FileInputStream("./mgampso-0.7548684662307973"));
        ArrayList<ArrayList<Pair<Integer, Integer>>> pairs = (ArrayList<ArrayList<Pair<Integer, Integer>>>) inputStream.readObject();
        System.out.println(pairs);
    }
}
