/**
 * Authors: Jesse Gerringer, Gleb Sklyr
 * Mails: gleb.sklyr@marquette.edu, Jesse.Gerringer@marquette.edu
 * Last edited: 12/11/2016 (mm/dd/yyyy) 7:24:PM
 * Description: Various methods that implement the Gibbs Sampler algorithm
 */
package gibbs_sampler;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.InputMismatchException;
import java.util.List;
import java.util.Random;
import java.util.Scanner;

public class Gibbs {
    
    /**
     * gets a random position based on the score
     * @param lWordsPr the scoe per word.  Must be in the words order
     * @return a randomposition index.  Higher probability to higher scores
     */
    protected static int getRandWordInd(List<Object> lWordsPr) {
        int word = randInt(0, lWordsPr.size()-1);
        double sum = 0;
        double[] frac = new double[lWordsPr.size()];
        double[] Norfrac = new double[lWordsPr.size()];
        for(Object score : lWordsPr) { // add up all scores
            sum = sum + (double)score;
        }
        int x = 0;
        for(Object score : lWordsPr) { // fill up the fraction array
            frac[x] = ((double)score/(double)sum)*100.0;
            x++;
        }
        for(int i = 0; i < frac.length; i++) { // normalize frac array
            Norfrac[i] = arraySum(frac, i);
        }
        int rand = Gibbs.randInt(0, 99); // random int
        for(int i = 0; i < Norfrac.length; i++){ // chrck of itsin range of i
            if(rand <= Norfrac[i]) {
                word = i;
                break; // found the random word!!
            }
        }
        // random word index (according to weight)
        return word;
    }
    
    /**
     * Sums an array till a certain point
     * @param array the array to sum
     * @param pos up to which point to sum
     * @return the sum of points in the array until pos
     */
    protected static double arraySum(double array[], int pos) {
        double sum = 0;
        for(int i = 0; i <= pos; i++) {
            sum = sum + array[i]; 
        }
        return sum;
    }
    
    /**
     * Prints the candidate motifs from each sequence in S to a file
     * @param outputStream the stream to the file
     */
    protected static void writeMotifToFile(PrintWriter outputStream) {
        for(int i = 0; i < Gibbs_Sampler.S.size(); i++) {
            for(int j = 0; j < Gibbs_Sampler.S.get(i).length(); j++) {
                if(Gibbs_Sampler.S.get(i).charAt(j) == 'A'||Gibbs_Sampler.S.get(i).charAt(j) == 'T'||Gibbs_Sampler.S.get(i).charAt(j) == 'C'||Gibbs_Sampler.S.get(i).charAt(j) == 'G') { // if its motif candidate
                    outputStream.println(Gibbs_Sampler.S.get(i).substring(j, (j+Gibbs_Sampler.l)));
                    break;
                }
            }
        }
    }
    
    /**
     * Prints the uppercased letters in each S set of sequences
     */
    protected static void printMotifOnly() {
        System.out.println("Candidate motif: ");
        for(int i = 0; i < Gibbs_Sampler.S.size(); i++) {
            for(int j = 0; j < Gibbs_Sampler.S.get(i).length(); j++) {
                if(Gibbs_Sampler.S.get(i).charAt(j) == 'A'||Gibbs_Sampler.S.get(i).charAt(j) == 'T'||Gibbs_Sampler.S.get(i).charAt(j) == 'C'||Gibbs_Sampler.S.get(i).charAt(j) == 'G') { // if its motif candidate
                    System.out.println(Gibbs_Sampler.S.get(i).substring(j, (j+Gibbs_Sampler.l)));
                    break;
                }
            }
        }
    }
    
    /**
     * Prints the concensus motif as from Theta
     */
    protected static void printConcensus() {
        System.out.println("Consensus Motif: ");
        for(int j = 0; j < Gibbs_Sampler.THETA_ATCG[0].length; j++) {
            String print = "A";
            double max = Gibbs_Sampler.THETA_ATCG[0][j];
            if(Gibbs_Sampler.THETA_ATCG[1][j] > max) {
                max = Gibbs_Sampler.THETA_ATCG[1][j];
                print = "T";
            } 
            if(Gibbs_Sampler.THETA_ATCG[2][j] > max) {
                max = Gibbs_Sampler.THETA_ATCG[2][j];
                print = "C";
            } 
            if(Gibbs_Sampler.THETA_ATCG[3][j] > max) {
                max = Gibbs_Sampler.THETA_ATCG[3][j];
                print = "G";
            }
            System.out.print(print);
        }
        System.out.println();
    }
    
    /**
     * Prints set S from Gibbs Sampler with the motif aligned
     */
    protected static void printSmotif() {
        System.out.println("S with aligned motif: ");
        int maxMotif = 0;
        List<String> motifPos = new ArrayList();
        // Find farthest motif
        for(String s : Gibbs_Sampler.S) { // for every string
            for(int i = 0; i < s.length(); i++) { // for every character
                if(s.charAt(i) == 'A'||s.charAt(i) == 'T'||s.charAt(i) == 'C'||s.charAt(i) == 'G') { // find first motif letter
                    motifPos.add("" + i); // add motif start pos
                    if(i > maxMotif) {
                    maxMotif = i;
                    }
                    break;
                }
            }
        }
        // Print with alignment
        int x = 0;
        for(String s : Gibbs_Sampler.S) {
            printSpaces(maxMotif-Integer.parseInt(motifPos.get(x)));
            System.out.print(s);
            System.out.println();
            x++;
        }
    }
    
    /**
     * prints some spaces in one line
     * @param numS the number of spaces to be printed
     */
    private static void printSpaces(int numS) {
        for(int i = 0; i < numS; i++) {
            System.out.print(" ");
        }
    }
    
    /**
     * Computes the sum of the PR scores for all set S
     * @return the sum of logs of the PR scores
     */
    protected static double logScoreSumS() {
        List<String> copy = Gibbs_Sampler.S;
        double sum = 0;
        for(String s : copy) {
            sum = sum + Math.log10(prScore(s));
        }
        return (double)sum;
    }
    
    /**
     * Computes the PR(Z|theta) / PR(Z|theta_zero)
     * @param x the string with uppercased motif candidate
     * @return the probability the motif is not in background quotient
     */
    protected static double prScore(String x) {
        int l = Gibbs_Sampler.l;
        int j = 0;
        double score = 1.0;
        double den = Gibbs_Sampler.THETA_0_ATCG[0]*Gibbs_Sampler.THETA_0_ATCG[1]*Gibbs_Sampler.THETA_0_ATCG[2]*Gibbs_Sampler.THETA_0_ATCG[3];
        while(j < l) { // for every motif letter
            for(int i = 0; i < x.length(); i++) { // for every letter in x
                if(x.charAt(i) == 'A'||x.charAt(i) == 'T'||x.charAt(i) == 'C'||x.charAt(i) == 'G') { // if its motif candidate
                    if(x.charAt(i+j) == 'A') {
                        score = score * Gibbs_Sampler.THETA_ATCG[0][j];
                        break; // next character
                    } else if(x.charAt(i+j) == 'T') {
                        score = score * Gibbs_Sampler.THETA_ATCG[1][j];
                        break; // next character
                    } else if(x.charAt(i+j) == 'C') {
                        score = score * Gibbs_Sampler.THETA_ATCG[2][j];
                        break; // next character
                    } else if(x.charAt(i+j) == 'G') {
                        score = score * Gibbs_Sampler.THETA_ATCG[3][j];
                        break; // next character
                    }
                }
            }
            j++;
        }
        return (double)((double)score/(double)den);
    }
    
    /**
     * prints the theta zero
     */
    protected static void printThetaZero(){System.out.println("Theta Zero: ");for(double i:Gibbs_Sampler.THETA_0_ATCG){System.out.print(i + ", ");}System.out.println();}
    
    /**
     * Generates the theta zero for the background
     */
    protected static void getThetaZero() {
        double[] copy = Gibbs_Sampler.THETA_0_ATCG;
        int count_a = 0;
        int count_t = 0;
        int count_c = 0;
        int count_g = 0;
        int count = 0;
        for(String s : Gibbs_Sampler.S) { // for every s
            for(int i = 0; i < s.length(); i++) { // for every char
                if(s.charAt(i) == 'a') {
                    count_a++;
                    count++;
                } else if(s.charAt(i) == 't') {
                    count_t++;
                    count++;
                } else if(s.charAt(i) == 'c') {
                    count_c++;
                    count++;
                } else if(s.charAt(i) == 'g') {
                    count_g++;
                    count++;
                }
            }
        }
        copy[0] = (double)count_a/(double)count;
        copy[1] = (double)count_t/(double)count;
        copy[2] = (double)count_c/(double)count;
        copy[3] = (double)count_g/(double)count;
        Gibbs_Sampler.THETA_0_ATCG = copy;
    }
    
    /**
     * print theta
     */
    protected static void printTheta(){System.out.println("Theta: ");for(int i=0;i<Gibbs_Sampler.THETA_ATCG.length;i++){for(int j=0;j<Gibbs_Sampler.THETA_ATCG[0].length;j++){System.out.print(" |"+Gibbs_Sampler.THETA_ATCG[i][j]);}System.out.println();}}
    
    /**
     * Loads the theta for the proposed motif in Gibbs Sampler Application
     */
    protected static void getTheta() {
        int l = Gibbs_Sampler.l;
        int countA = 0;
        int countT = 0;
        int countC = 0;
        int countG = 0;
        double[][] thetaATCG = new double[4][l];
        int j = 0;
        while(j < l) // for every motif char
        {
        for(String s : Gibbs_Sampler.S) { // for every s
            for(int i = 0; i < s.length(); i++) { // for every letter
                if(s.charAt(i) == 'A'||s.charAt(i) == 'T'||s.charAt(i) == 'C'||s.charAt(i) == 'G')
                {
                if(s.charAt(i + j) == 'A') {
                    countA++;
                    break;
                } else if(s.charAt(i + j) == 'T') {
                    countT++;
                    break;
                } else if(s.charAt(i + j) == 'C') {
                    countC++;
                    break;
                } else if(s.charAt(i + j) == 'G') {
                    countG++;
                    break;
                }
                }
            }
        }
        thetaATCG[0][j] = (double)countA/Gibbs_Sampler.S.size();
        countA = 0;
        thetaATCG[1][j] = (double)countT/Gibbs_Sampler.S.size();
        countT = 0;
        thetaATCG[2][j] = (double)countC/Gibbs_Sampler.S.size();
        countC = 0;
        thetaATCG[3][j] = (double)countG/Gibbs_Sampler.S.size();
        countG = 0;
        j++;
        }
        Gibbs_Sampler.THETA_ATCG = thetaATCG;
    }
    
    /**
     * print S from Gibbs Sampler for debug purposes
     */
    protected static void printS(){System.out.println("S: ");Gibbs_Sampler.S.forEach((s)->{System.out.println(s);});}
    
    /**
     * uppercases a random Z big (set of motifs) in S set of string that
     * is obtained in Gibbs_Sampler application
     */
    protected static void uppercaseRnadZ() {
        int l = Gibbs_Sampler.l;
        List<String> copy = Gibbs_Sampler.S;
        List<String> Z = new ArrayList();
        for(String s : copy) {
            int randStart = getRandMotifPos(s); // randon z little statrt pos
            String begin = s.substring(0, randStart);
            String motif = s.substring(randStart, randStart + l).toUpperCase();
            String rest = "";
            if(randStart + l != s.length()) { // if the motif is not the end ofo s
                rest = s.substring(randStart + l, s.length());
            }
            Z.add(begin + motif + rest);
        }
        Gibbs_Sampler.S = Z;
    }
    
    /**
     * get a random start position for a length l motif
     * l is obtained from Gibbs Sampler application
     * @param S the string to generate start position for a motif
     * @return the random motif start position
     */
    private static int getRandMotifPos(String S) {
        int motifL = Gibbs_Sampler.l;
        int maxStart = S.length() - motifL; // maximum start position
        int startPos = randInt(0, maxStart);
        return startPos;
    }
    
    /**
     * Generates a random number in a specified range (including both extremes)
     * @param min the min value desired (including)
     * @param max the max value desired (including)
     * @return 
     */
    protected static int randInt(int min, int max) {
        Random random = new Random();
        int rand;
        rand = random.nextInt(max - min + 1) + min;
        return rand;
    }
    
    /**
    * Interacts with the user until the user enters a positive non zero integer
    * The method will check for malicious input and display appropriate
    * errors for the user to retry. Loads the integer to Gibbs Sampler application
    * @param prompt What should we ask the user
    */
    protected static void acceptIntNotZero(String prompt) {
        int user = 0;
        Scanner keyboard = new Scanner(System.in);
        while(true) { // keep asking the user
            user = 0;
            try { // check for malicious input
                System.out.print(prompt);
                user = keyboard.nextInt();
            } catch (InputMismatchException ex) {
                System.out.println("Must enter an integer! Try again...");
                keyboard.next();
                continue; // re-loop
            }
            if(user <= 0) // input is good but negative or zero
                System.out.println("A value greater than zero is required");
            else {
                break; // valid input entered
            }
        }
        Gibbs_Sampler.l = user;
    }
    
    /**
     * Loads a file specified by the user into the S set of strings
     * in the Gibbs_Sampler application.  Lower cases everything!
     * @param prompt the specific file prompt
     */
    protected static void getS(String prompt) {
        List<String> user = new ArrayList();
        Scanner inputStream = promptFile(prompt);
        while(inputStream.hasNextLine()) {
            user.add(inputStream.nextLine().trim().toLowerCase());
        }
        Gibbs_Sampler.S = user;
    }
    
    /**
    * upon user input this method tries to find an existing file
    * to write to.  If fails, prompts again until existing file
    * name received
    * @param prompt the specific file asked for
    * @return Scanner object ready to scan from existing file
    */
    private static Scanner promptFile(String prompt) {
        //promt user for existing file only
        while(true) {
            try {
                System.out.print(prompt);
                Scanner console = new Scanner(System.in);
                String inputFileName = console.next();
                return new Scanner(new FileInputStream(inputFileName));
            }
            catch (FileNotFoundException ex) {
                System.out.println("File does not exist.  Enter existing file.  Error "+ ex.getLocalizedMessage());
            }
        }
    }
    
    /**
    * prompts for output file name from user.  If file name is in
    * use, prompts again until unused name received
    * @return the file to write to
    */
    public static File promptOutFile() {
        while(true) {
            Scanner console = new Scanner(System.in);
            System.out.println("enter output file name for motif(out.txt)");
            String outputFileName = console.next();
            if(new File(outputFileName).exists())
                System.out.println("File name exits.  Try a different one");
            else
                return new File(outputFileName);
        }
    }
    
    /**
     * tries to connect a PrintWriter object to an output
     * file by searching by name and if file with such
     * name does not exist - throws an exception
     * @param outputFileName String name of the file to open
     * @return a PrintWriter object ready to write to file
     */
    public static PrintWriter outputStreamer(String outputFileName) {
        while(true) {
            try  {
                return new PrintWriter(outputFileName);
            }
            catch (FileNotFoundException ex) {
                System.out.println("no valid output file name");
            }
        }
    } 
}