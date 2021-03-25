package util;

import org.kamranzafar.commons.cloner.ObjectCloner;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

public class Helper {
    /**
     * getColLength
     * Returns an integer 1.2 times greater than the longest sequence
     * in the given sequences.
     *
     * @param seqs MSA util.Sequences
     * @return column length
     */
    public static int getColLength(String[] seqs) {
        int length = 0;
        for (String seq : seqs) {
            if (seq.length() > length)
                length = seq.length();
        }
        return (int) Math.ceil(length * 1.2);
    }

    /**
     * An aggregated function.
     *
     * @param bitmatrix
     * @param seq
     * @param w1
     * @param w2
     * @param checkInfeasibility
     * @param ops
     * @return
     */
    public static double aggregatedFunction(int[][] bitmatrix,
                                            String[] seq,
                                            double w1,
                                            double w2,
                                            boolean checkInfeasibility,
                                            Operator[] ops) {
        if (checkInfeasibility && ops == null)
            throw new IllegalArgumentException("Checking for infeasibility but ops are null!");

        if (checkInfeasibility && infeasible(bitmatrix, seq, ops))
            return Double.MIN_VALUE;

        String[] strings = bitsToStrings(bitmatrix, seq);

        int nMax = maxNumOfIndels(bitmatrix, seq); // total number of indels

        int nI = numOfInsertedIndels(bitmatrix, seq); // number of indels before last char

        return (w1 * numOfAlignedChars(strings)) + (w2 * (nMax - nI));
    }

    /**
     * Counts the number of indels that are found
     * before the last character in a sequence.
     */
    public static int numOfInsertedIndels(int[][] bitmatrix,
                                          String[] seq) {
        /*
        Remember: bit 0 means character from sequence
                     bit 1 means indel
         */

        int count = 0;

        for (int i = 0; i < seq.length; i++) { // Go through each bitstring in bitmatrix
            int tmp = 0;
            boolean hitLastChar = false;
            for (int bit : bitmatrix[i]) { // Go through each bit in bitstring
                if (bit == 0) { // we hit a character
                    tmp += 1;
                    if (tmp == seq[i].length())
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
    public static int maxNumOfIndels(int[][] bitmatrix,
                                     String[] seq) {
        int count = 0;
        for (int i = 0; i < bitmatrix.length; i++)
            count += bitmatrix[i].length - seq[i].length();
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
    public static boolean infeasible(int[][] bitmatrix,
                                     String[] seq,
                                     Operator[] ops) {
        for (int i = 0; i < seq.length; i++) {
            int count = 0;
            for (int bit : bitmatrix[i])
                if (bit == 0)
                    count += 1;

            for (Operator op : ops)
                if (op.operation(count, seq[i].length()))
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
    public static String[] bitsToStrings(int[][] bitmatrix,
                                         String[] seq) {
        String[] result = new String[seq.length];
        int i = 0;
        for (int[] bitlist : bitmatrix) {
            int j = 0;
            StringBuilder tmp = new StringBuilder();
            for (int bit : bitlist) {
                if (bit == 0 && j < seq[i].length()) {
                    tmp.append(seq[i].charAt(j));
                    j++;
                } else
                    tmp.append("-");
            }
            result[i] = tmp.toString();
            i++;
        }
        return result;
    }

    /**
     * Counts the number of aligned characters in a list of strings.
     * <p>
     * Important: Assumes that each string is of the same length.
     */
    public static int numOfAlignedChars(String[] strings) {
        if (strings.length < 1)
            throw new IllegalArgumentException("List of strings is empty");
        else if (strings.length == 1) {
            System.out.println("Warning: only 1 string in numOfAlignedChars function!");
            return 0;
        }

        int result = 0;
        HashMap<Character, Integer> charList = new HashMap<>();

        /*
            The charList is a mapping, of each character to
            the amount of matchings that character has.
         */
        for (int i = 0; i < strings[0].length(); i++) {
            HashMap<Character, Integer> tmpMapping = new HashMap<>(); // Temp before adding to main charList
            // Fill out the tmpMapping for all characters in column i
            for (String str : strings) {
                Character c = str.charAt(i);
                if (c != '-') {
                    if (tmpMapping.containsKey(c))
                        tmpMapping.put(c, tmpMapping.get(c) + 1);
                    else
                        tmpMapping.put(c, 1);
                }
            }

            // Now add the tmpMapping into main charList
            for (Character c : tmpMapping.keySet()) {
                // if there's more than one of that character in the col, there's a matching
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
    public static int[][] genBitMatrix(double[] pos,
                                       String[] seq,
                                       int colLength) {
        int[][] bitmatrix = new int[seq.length][colLength];

        for (int i = 0; i < seq.length; i++) {
            for (int j = 0; j < colLength; j++) {

                double value = gen(((i * colLength) + j) * 0.1, pos[0], pos[1], pos[2], pos[3]);
                if (value > 0)
                    bitmatrix[i][j] = 1;
                else
                    bitmatrix[i][j] = 0;
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

    /**
     * Checks if one bitmatrix dominates another.
     * <p>
     * This means that for each FitnessFunction,
     * a bitmatrix must be at least as good as
     * the other in all functions, and in at least
     * one function it must be better.
     *
     * @param seq
     * @param bm1
     * @param bm2
     * @return
     */
    private static boolean dominates(String[] seq,
                                     int[][] bm1,
                                     int[][] bm2,
                                     FitnessFunction[] f) {
        boolean betterInAtLeastOne = false;

        for (FitnessFunction ff : f) {
            double res1 = ff.calculate(bm1, seq);
            double res2 = ff.calculate(bm2, seq);

            if (res1 <= res2)
                if (res1 < res2)
                    betterInAtLeastOne = true;
                else
                    return false;
        }

        return betterInAtLeastOne;
    }

    /**
     * Checks that two bitmatrixes are the same.
     *
     * @param x
     * @param y
     * @return
     */
    private static boolean theSame(int[][] x,
                                   int[][] y) {
        for (int i = 0; i < x.length; i++) {
            if (x[i].length != y[i].length)
                return false;

            for (int j = 0; j < x[i].length; j++) {
                if (x[i][j] != y[i][j])
                    return false;
            }
        }

        return true;
    }

    /**
     * Solution x is attempted to be added into the
     * archive if it is not dominated.
     * <p>
     * If after addition to the archive, and previous
     * solutions are now dominated by this new addition,
     * then those dominated solutions are removed.
     * <p>
     * Also, if the archive is full then the most crowded
     * solution is removed.
     *
     * @param seq
     * @param sArchive
     * @param x
     * @param f
     * @return
     */
    public static void addToArchiveB(String[] seq,
                                     ArrayList<Pair<int[][], Double>> sArchive,
                                     int[][] x,
                                     FitnessFunction[] f,
                                     int limit) {
        ArrayList<Pair<int[][], Double>> aDominated = new ArrayList<>();

        for (Pair<int[][], Double> s : sArchive) {
            // archive entry dominates x
            if (dominates(seq, s.getFirst(), x, f))
                return;
            else if (dominates(seq, x, s.getFirst(), f)) // x dominates archive entry
                aDominated.add(s);

            // x and archive entry are the exact same
            if (theSame(x, s.getFirst()))
                return;
        }

        /*
         When function enters here, x is not dominated by any
         entry in the archive, but may have dominated some
         entries. So add it to archive, but remove all the
         previous entries that are dominated by it.

         If the limit of the archive has been reached, then
         calculate crowding distances and remove the most
         crowded solution.
         */

        sArchive.removeAll(aDominated);

        /*
             If this is true then the archive is completely full so
             remove most crowded solution.
         */
        if (sArchive.size() == limit) {
            calculateCrowdingDistancesB(seq, sArchive, f);

            Pair<int[][], Double> mostCrowdedEntry = sArchive.get(0);

            for (Pair<int[][], Double> entry : sArchive) {
                if (entry.getSecond() < mostCrowdedEntry.getSecond())
                    mostCrowdedEntry = entry;
            }

            sArchive.remove(mostCrowdedEntry);
        }

        ObjectCloner<int[][]> cloner = new ObjectCloner<>();
        sArchive.add(new Pair<>(cloner.deepClone(x), 0.0));
    }

    public static void addToArchiveA(String[] seq,
                                     ArrayList<Triplet<double[], int[][], Double>> sArchive,
                                     Pair<double[], int[][]> x,
                                     FitnessFunction[] f,
                                     int limit) {
        ArrayList<Triplet<double[], int[][], Double>> aDominated = new ArrayList<>();

        for (Triplet<double[], int[][], Double> s : sArchive) {
            // archive entry dominates x
            if (dominates(seq, s.getSecond(), x.getSecond(), f))
                return;
            else if (dominates(seq, x.getSecond(), s.getSecond(), f)) // x dominates archive entry
                aDominated.add(s);

            // x and archive entry are the exact same
            if (theSame(x.getSecond(), s.getSecond()))
                return;
        }

        /*
         When function enters here, x is not dominated by any
         entry in the archive, but may have dominated some
         entries. So add it to archive, but remove all the
         previous entries that are dominated by it.

         If the limit of the archive has been reached, then
         calculate crowding distances and remove the most
         crowded solution.
         */

        sArchive.removeAll(aDominated);

        /*
             If this is true then the archive is completely full so
             remove most crowded solution.
         */
        if (sArchive.size() == limit) {
            calculateCrowdingDistancesA(seq, sArchive, f);

            Triplet<double[], int[][], Double> mostCrowdedEntry = sArchive.get(0);

            for (Triplet<double[], int[][], Double> entry : sArchive) {
                if (entry.getThird() < mostCrowdedEntry.getThird())
                    mostCrowdedEntry = entry;
            }

            sArchive.remove(mostCrowdedEntry);
        }

        ObjectCloner<int[][]> cloner = new ObjectCloner<>();
        ObjectCloner<double[]> cloner2 = new ObjectCloner<>();
        sArchive.add(new Triplet<>(cloner2.deepClone(x.getFirst()), cloner.deepClone(x.getSecond()), 0.0));
    }

    /**
     * Get a guide from the archive.
     *
     * @param seq      sequences
     * @param sArchive the archive
     * @param f        list of functions
     * @param k        number of elements for tournament selection
     * @return archive value as guide
     */
    public static int[][] archiveGuideB(String[] seq,
                                        ArrayList<Pair<int[][], Double>> sArchive,
                                        FitnessFunction[] f,
                                        int k) {

        calculateCrowdingDistancesB(seq, sArchive, f);

        int index = ThreadLocalRandom.current().nextInt(0, sArchive.size());

        for (int i = 0; i < k; i++) {
            int tmp = ThreadLocalRandom.current().nextInt(0, sArchive.size());
            if (sArchive.get(index).getSecond() > sArchive.get(tmp).getSecond())
                index = tmp;
        }

        return sArchive.get(index).getFirst();
    }

    public static double[] archiveGuideA(String[] seq,
                                         ArrayList<Triplet<double[], int[][], Double>> sArchive,
                                         FitnessFunction[] f,
                                         int k) throws Exception {

        if(sArchive.size() == 0)
            return null;

        calculateCrowdingDistancesA(seq, sArchive, f);

        int index = ThreadLocalRandom.current().nextInt(0, sArchive.size());

        for (int i = 0; i < k; i++) {
            int tmp = ThreadLocalRandom.current().nextInt(0, sArchive.size());
            if (sArchive.get(index).getThird() > sArchive.get(tmp).getThird())
                index = tmp;
        }

        return sArchive.get(index).getFirst();
    }

    /**
     * Crowding distance calculation for archive
     * entry removal.
     * <p>
     * Uses the updated crowding distance calculation from:
     *
     * @param f
     */
    private static ArrayList<Pair<int[][], Double>> calculateCrowdingDistancesB(String[] seq,
                                                                                ArrayList<Pair<int[][], Double>> sArchive,
                                                                                FitnessFunction[] f) {

        for (FitnessFunction ff : f) {
            // Sort sArchive in ascending order
            Comparator<Pair<int[][], Double>> c = (o1, o2) -> ff.calculate(o1.getFirst(), seq) > ff.calculate(o2.getFirst(), seq) ? 1 : 0;
            sArchive.sort(c);

            // Set the boundary distances to infinite (max value)
            sArchive.get(0).setSecond(Double.MAX_VALUE);
            sArchive.get(sArchive.size() - 1).setSecond(Double.MAX_VALUE);

            for (int i = 1; i < sArchive.size() - 1; i++) {
                double curr = sArchive.get(i).getSecond();
                double next = ff.calculate(sArchive.get(i + 1).getFirst(), seq);
                double prev = ff.calculate(sArchive.get(i - 1).getFirst(), seq);
                sArchive.get(i).setSecond(curr + next - prev);
            }
        }

        // This will sort uniqueDistancedFitnesses based on distance (ascending)
        Comparator<Pair<int[][], Double>> c = (o1, o2) -> o1.getSecond() > o2.getSecond() ? 1 : 0;
        sArchive.sort(c);

        return sArchive;
    }

    private static ArrayList<Triplet<double[], int[][], Double>> calculateCrowdingDistancesA(String[] seq,
                                                                                             ArrayList<Triplet<double[], int[][], Double>> sArchive,
                                                                                             FitnessFunction[] f) {

        for (FitnessFunction ff : f) {
            // Sort sArchive in ascending order
            Comparator<Triplet<double[], int[][], Double>> c = (o1, o2) -> ff.calculate(o1.getSecond(), seq) > ff.calculate(o2.getSecond(), seq) ? 1 : 0;
            sArchive.sort(c);

            // Set the boundary distances to infinite (max value)
            sArchive.get(0).setThird(Double.MAX_VALUE);
            sArchive.get(sArchive.size() - 1).setThird(Double.MAX_VALUE);

            for (int i = 1; i < sArchive.size() - 1; i++) {
                double curr = sArchive.get(i).getThird();
                double next = ff.calculate(sArchive.get(i + 1).getSecond(), seq);
                double prev = ff.calculate(sArchive.get(i - 1).getSecond(), seq);
                sArchive.get(i).setThird(curr + next - prev);
            }
        }

        // This will sort uniqueDistancedFitnesses based on distance (ascending)
        Comparator<Triplet<double[], int[][], Double>> c = (o1, o2) -> o1.getThird() > o2.getThird() ? 1 : 0;
        sArchive.sort(c);

        return sArchive;
    }

    public static double average(double[] arr) {
        double total = 0;
        for (double d : arr)
            total += d;
        return total / arr.length;
    }

    public static double average(int[] arr) {
        double total = 0;
        for (double d : arr)
            total += d;
        return total / arr.length;
    }

    public static double stdev(double[] arr) {
        double mean = average(arr);
        double sum = 0;
        for (double d : arr) {
            sum += Math.pow(d - mean, 2);
        }
        return Math.sqrt(sum / arr.length);
    }

    public static double stdev(int[] arr) {
        double mean = average(arr);
        double sum = 0;
        for (double d : arr) {
            sum += Math.pow(d - mean, 2);
        }
        return Math.sqrt(sum / arr.length);
    }

    public static void printInfo(String title, double highest, double lowest, double[] list) {
        System.out.println("Best " + title + ": " + highest);
        System.out.println("Worst " + title + ": " + lowest);
        System.out.println("Average " + title + ": " + Helper.average(list));
        System.out.println("St. Dev. " + title + ": " + Helper.stdev(list));
        System.out.println();
    }

    public static void printInfo(String title, int highest, int lowest, int[] list) {
        System.out.println("Highest " + title + ": " + highest);
        System.out.println("Lowest " + title + ": " + lowest);
        System.out.println("Average " + title + ": " + Helper.average(list));
        System.out.println("St. Dev. " + title + ": " + Helper.stdev(list));
        System.out.println();
    }

    public static void printInfo(String title, double[] highest, double[] lowest, double[][] list) {
        for (int i = 0; i < highest.length; i++) {
            System.out.println("Highest " + title + "[" + i + "]: " + highest[i]);
            System.out.println("Lowest " + title + "[" + i + "]: " + lowest[i]);
            System.out.println("Average " + title + "[" + i + "]: " + Helper.average(list[i]));
            System.out.println("St. Dev. " + title + "[" + i + "]: " + Helper.stdev(list[i]));
            System.out.println();
        }
    }

    public static void printInfo(String title, int[] highest, int[] lowest, int[][] list) {
        for (int i = 0; i < highest.length; i++) {
            System.out.println("Highest " + title + "[" + i + "]: " + highest[i]);
            System.out.println("Lowest " + title + "[" + i + "]: " + lowest[i]);
            System.out.println("Average " + title + "[" + i + "]: " + Helper.average(list[i]));
            System.out.println("St. Dev. " + title + "[" + i + "]: " + Helper.stdev(list[i]));
            System.out.println();
        }
    }
}
