package util;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class TFAParser {
    public static String[] seq1;
    public static String[] seq2;
    public static String[] seq3;
    public static String[] seq4;
    public static String[] seq5;
    public static String[] seq6;
    public static String[] seq7;
    public static String[] seq8;
    public static String[] seq9;
    public static String[] seq10;
    public static String[] seq11;
    public static String[] seq12;
    public static String[] seq13;

    static {
        try {
            seq1 = getStrings(new File("./bb3_release/RV11/BB11001.tfa"));
            seq2 = getStrings(new File("./bb3_release/RV11/BB11002.tfa"));
            seq3 = getStrings(new File("./bb3_release/RV11/BB11005.tfa"));
            seq4 = getStrings(new File("./bb3_release/RV11/BB11009.tfa"));
            seq5 = getStrings(new File("./bb3_release/RV11/BB11013.tfa"));
            seq6 = getStrings(new File("./bb3_release/RV11/BB11026.tfa"));
            seq7 = getStrings(new File("./bb3_release/RV12/BB12003.tfa"));
            seq8 = getStrings(new File("./bb3_release/RV12/BB12005.tfa"));
            seq9 = getStrings(new File("./bb3_release/RV12/BB12006.tfa"));
            seq10 = getStrings(new File("./bb3_release/RV12/BB12009.tfa"));
            seq11 = getStrings(new File("./bb3_release/RV12/BB12014.tfa"));
            seq12 = getStrings(new File("./bb3_release/RV12/BB12020.tfa"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static String[] getStrings(File file) throws FileNotFoundException {
        Scanner input = new Scanner(file);
        ArrayList<String> strings = new ArrayList<>();
        StringBuilder newString = new StringBuilder();
        while (input.hasNextLine()) {
            String line = input.nextLine();
            if (line.charAt(0) == '>') {
                if (!newString.toString().isEmpty())
                    strings.add(newString.toString());
                newString = new StringBuilder();
            } else {
                newString.append(line);
            }
        }

        if(!newString.toString().isEmpty())
            strings.add(newString.toString());

        String[] newStrings = new String[strings.size()];
        for (int i = 0; i < newStrings.length; i++)
            newStrings[i] = strings.get(i);

        return newStrings;
    }

    public static void main(String[] args) throws FileNotFoundException {
        // RV11 Set
        File file = new File("./bb3_release/RV11/BB11001.tfa");
        File file2 = new File("./bb3_release/RV11/BB11002.tfa");
        File file3 = new File("./bb3_release/RV11/BB11005.tfa");
        File file4 = new File("./bb3_release/RV11/BB11009.tfa");
        File file5 = new File("./bb3_release/RV11/BB11013.tfa");
        File file6 = new File("./bb3_release/RV11/BB11026.tfa");
        System.out.println(Arrays.toString(getStrings(file)));
        System.out.println(Arrays.toString(getStrings(file2)));
        System.out.println(Arrays.toString(getStrings(file3)));
        System.out.println(Arrays.toString(getStrings(file4)));
        System.out.println(Arrays.toString(getStrings(file5)));
        System.out.println(Arrays.toString(getStrings(file6)));

        // RV12 Set
        File file7 = new File("./bb3_release/RV12/BB12003.tfa");
        File file8 = new File("./bb3_release/RV12/BB12005.tfa");
        File file9 = new File("./bb3_release/RV12/BB12006.tfa");
        File file10 = new File("./bb3_release/RV12/BB12009.tfa");
        File file11 = new File("./bb3_release/RV12/BB12014.tfa");
        File file12 = new File("./bb3_release/RV12/BB12020.tfa");
        System.out.println(Arrays.toString(getStrings(file7)));
        System.out.println(Arrays.toString(getStrings(file8)));
        System.out.println(Arrays.toString(getStrings(file9)));
        System.out.println(Arrays.toString(getStrings(file10)));
        System.out.println(Arrays.toString(getStrings(file11)));
        System.out.println(Arrays.toString(getStrings(file12)));

        String[][] strs = {seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8, seq9, seq10, seq11, seq12};

        for(String[] s : strs) {
            System.out.print(s.length + " ");
            for(String ss : s)
                System.out.print(ss.length() + " ");
            System.out.println();
        }
    }
}
