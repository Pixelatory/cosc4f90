package util;

import org.kamranzafar.commons.cloner.ObjectCloner;

import java.util.*;

public class Helper {

    // Standard operators for use in algorithms
    public static Operator lt = (a, b) -> a < b;
    public static Operator gt = (a, b) -> a > b;
    public static Operator eq = (a, b) -> a == b;

    /**
     * getColLength
     * Returns an integer 1.2 times greater than the longest sequence
     * in the given sequences.
     *
     * @param seqs MSA util.Sequences
     * @return column length
     */
    public static int getColLength(ArrayList<String> seqs) {
        int length = 0;
        for (String seq : seqs) {
            if (seq.length() > length)
                length = seq.length();
        }
        return (int) Math.ceil(length * 1.2);
    }

    public static <E> ArrayList<ArrayList<E>> copyArray(ArrayList<ArrayList<E>> arr) {
        ArrayList<ArrayList<E>> tmp = new ArrayList<>();
        ObjectCloner<ArrayList<E>> cloner = new ObjectCloner<>();
        for (ArrayList<E> subArr : arr) {
            tmp.add(cloner.deepClone(subArr));
        }
        return tmp;
    }

    public static double aggregatedFunction(ArrayList<ArrayList<Integer>> bitmatrix,
                                            ArrayList<String> seq,
                                            double w1,
                                            double w2,
                                            boolean checkInfeasibility,
                                            ArrayList<Operator> ops) {
        if (checkInfeasibility && ops == null)
            throw new IllegalArgumentException("Checking for infeasibility but ops are null!");

        if (checkInfeasibility && infeasible(bitmatrix, seq, ops))
            return Double.MIN_VALUE;

        ArrayList<String> strings = bitsToStrings(bitmatrix, seq);

        int nMax = maxNumOfIndels(bitmatrix, seq); // total number of indels

        int nI = numOfInsertedIndels(bitmatrix, seq); // number of indels before last char

        return (w1 * numOfAlignedChars(strings)) + (w2 * (nMax - nI));
    }

    /**
     * Counts the number of indels that are found before the last character in a sequence.
     */
    public static int numOfInsertedIndels(ArrayList<ArrayList<Integer>> bitmatrix,
                                          ArrayList<String> seq) {
        /* Remember: bit 0 means character from sequence
                     bit 1 means indel
         */

        int count = 0;

        for (int i = 0; i < seq.size(); i++) { // Go through each bitstring in bitmatrix
            int tmp = 0;
            boolean hitLastChar = false;
            for (int bit : bitmatrix.get(i)) { // Go through each bit in bitstring
                if (bit == 0) { // we hit a character
                    tmp += 1;
                    if (tmp == seq.get(i).length())
                        hitLastChar = true;
                } else // we hit an indel
                    count += 1;

                if (hitLastChar)
                    break;
            }
        }

        return count;
    }

    /**
     * Counts the total number of indels according to the difference between
     * bitmatrix row length and sequence length.
     * <p>
     * Does not consider individual bits.
     */
    public static int maxNumOfIndels(ArrayList<ArrayList<Integer>> bitmatrix,
                                     ArrayList<String> seq) {
        int count = 0;
        for (int i = 0; i < bitmatrix.size(); i++)
            count += bitmatrix.get(i).size() - seq.get(i).length();
        return count;
    }

    /**
     * Determines if a bitmatrix and sequence combination is infeasible.
     * <p>
     * Counts the number of 0 bits in each row of the bitmatrix, and if it
     * doesn't match up with an operator comparison between the 0 bit count
     * and the length of the sequence, returns True for infeasibility.
     * <p>
     * The operator can be anything within the parameter format, however
     * usually it is >, <, =. For easy access to these use the operator
     * module.
     * <p>
     * Ex. Say ops is [>, <], then bitmatrix is infeasible when # of 0's > # of chars,
     * and when # of 0's < # of chars. Thus it will only be feasible when # of 0's = # of chars.
     */
    public static boolean infeasible(ArrayList<ArrayList<Integer>> bitmatrix,
                                     ArrayList<String> seq,
                                     ArrayList<Operator> ops) {
        for (int i = 0; i < seq.size(); i++) {
            int count = 0;
            for (int bit : bitmatrix.get(i))
                if (bit == 0)
                    count += 1;

            for (Operator op : ops)
                if (op.operation(count, seq.get(i).length()))
                    return true;
        }

        return false;
    }

    /**
     * Converts a list of sequences into a list of strings with indels, according to the bitmatrix provided.
     * <p>
     * Note: A bit of 0 means a character, and a bit of 1 means indel.
     * However, if the number of 0 bits exceeds the length of the sequence,
     * indels will be inserted instead.
     */
    public static ArrayList<String> bitsToStrings(ArrayList<ArrayList<Integer>> bitmatrix,
                                                  ArrayList<String> seq) {
        ArrayList<String> result = new ArrayList<>();
        int i = 0;
        for (ArrayList<Integer> bitlist : bitmatrix) {
            int j = 0;
            StringBuilder tmp = new StringBuilder();
            for (int bit : bitlist) {
                if (bit == 0 && j < seq.get(i).length()) {
                    tmp.append(seq.get(i).charAt(j));
                    j += 1;
                } else
                    tmp.append("-");
            }
            result.add(tmp.toString());
            i += 1;
        }
        return result;
    }

    /**
     * Counts the number of aligned characters in a list of strings.
     * <p>
     * Important: Assumes that each string is of the same length.
     */
    public static int numOfAlignedChars(ArrayList<String> strings) {
        if (strings.size() < 1)
            throw new IllegalArgumentException("List of strings is empty");
        else if (strings.size() == 1) {
            System.out.println("Warning: only 1 string in numOfAlignedChars function!");
            return 0;
        }

        int result = 0;
        HashMap<Character, Integer> charList = new HashMap<>();

        /* The charList is a mapping, of each character to
         */
        for (int i = 0; i < strings.get(0).length(); i++) {
            HashMap<Character, Integer> tmpMapping = new HashMap<>(); // Temp before adding to main charList
            for (String str : strings) { // Fill out the tmpMapping for all characters in column i
                Character c = str.charAt(i);
                if (c != '-') {
                    if (tmpMapping.containsKey(c))
                        tmpMapping.put(c, tmpMapping.get(c) + 1);
                    else
                        tmpMapping.put(c, 1);
                }
            }

            for (Character c : tmpMapping.keySet()) { // Now add the tmpMapping into main charList
                if (tmpMapping.get(c) > 1) {
                    if (charList.containsKey(c))
                        charList.put(c, charList.get(c) + tmpMapping.get(c));
                    else
                        charList.put(c, tmpMapping.get(c));
                }
            }
        }

        for (Integer i : charList.values())
            result += i;

        return result;
    }

    /**
     * The classic sigmoid function.
     */
    public static double Sigmoid(double x) {
        return 1 / (1 + Math.exp(-x));
    }

    /**
     * Generating a bitMatrix from the
     * angular modulated generation function.
     * <p>
     * Values are sampled from 0, 1, ..., n-1
     *
     * @param pos
     * @param seq
     * @param colLength
     * @return
     */
    public static ArrayList<ArrayList<Integer>> genBitMatrix(ArrayList<Double> pos,
                                                             ArrayList<String> seq,
                                                             int colLength) {
        ArrayList<ArrayList<Integer>> bitmatrix = new ArrayList<>();

        for (int i = 0; i < seq.size(); i++) {
            bitmatrix.add(new ArrayList<>());
            for (int j = 0; j < colLength; j++) {
                double value = gen((i * colLength) + j, pos.get(0), pos.get(1), pos.get(2), pos.get(3));
                if (value > 0)
                    bitmatrix.get(i).add(1);
                else
                    bitmatrix.get(i).add(0);
            }
        }

        return bitmatrix;
    }

    /**
     * Angular Modulated Generation Function
     *
     * @param x
     * @param a
     * @param b
     * @param c
     * @param d
     * @return float
     */
    private static double gen(double x,
                       double a,
                       double b,
                       double c,
                       double d) {
        return Math.sin(2 * Math.PI * (x - a) * b * Math.cos(2 * Math.PI * c * (x - a))) + d;
    }
}
