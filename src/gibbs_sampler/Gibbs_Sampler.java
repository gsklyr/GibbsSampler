/**
 * Authors: Jesse Gerringer, Gleb Sklyr
 * Mails: gleb.sklyr@marquette.edu, Jesse.Gerringer@marquette.edu
 * Last edited: Last edited: 12/11/2016 (mm/dd/yyyy) 7:24:PM
 * Description: This program reads a text file containing sequences with possible
 * motif and a length l from the user.  It runs the Gibbs Sampler algorithm
 * on the sequence to find a candidate motif.  It then prints out the candidate
 * motif sequences from each sequence to a text file defined by the user.  There
 * are various options for console printouts in the end of main.
 */
package gibbs_sampler;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class Gibbs_Sampler {

    protected static List<String> S = new ArrayList(); // Set of sequences
    protected static int l; // Length of z - little (one motif)
    protected static double[][] THETA_ATCG;
    protected static double[] THETA_0_ATCG = new double[4];
    
    public static void main(String[] args) {
        File out;
        
        // Obtain Gibbs Sampler parameters S and l (10.9 ensure + require)
        Gibbs.getS("Please enter the file path for file containing S: ");
        Gibbs.acceptIntNotZero("Please enter the l of desired motif: ");
        
        // Generate Z in S of length l for each z (10.9 1)
        Gibbs.uppercaseRnadZ();
        
        // Get theta for Z and theta zero for background (10.9 2)
        Gibbs.getTheta();
        Gibbs.getThetaZero();
        
        // Get initial motif score 
        double score = Gibbs.logScoreSumS();
        double score2;
        
        int counter = 0; // counts the number of times the log score does not change
        // Sampler loop (10.9 3)
        while(true) {
            int randSeq = Gibbs.randInt(0, S.size()-1); // select random seq (10.9 4)
            S.set(randSeq, S.get(randSeq).toLowerCase()); // delete the word (10.9 4)
            Gibbs.getTheta(); // (10.9 5)
            Gibbs.getThetaZero(); // (10.9 5)
            List<Object> lWordsPr = new ArrayList(); // Pr score for all l words in Si
            List<Object> lWords = new ArrayList(); // all l words in Si
            // Figure 10.9 (6)
            for(int i = 0; i < S.get(randSeq).length() - l + 1; i++) { // for every l word in Si
                String begin = S.get(randSeq).substring(0, i);
                String motif = S.get(randSeq).substring(i, i + l).toUpperCase();
                String rest = "";
                if(i + l != S.get(randSeq).length()) { // if the motif is not the end of s
                    rest = S.get(randSeq).substring(i + l, S.get(randSeq).length());
                }
                S.set(randSeq, (begin + motif + rest)); // set a new word
                double newScore = Gibbs.prScore(S.get(randSeq));
                lWordsPr.add(newScore); // add the words score
                lWords.add(S.get(randSeq)); // add the word
                // delete the word
                S.set(randSeq, S.get(randSeq).toLowerCase());
            }
            
            // Choose the word of length l randomly depending on weight (10.9 7)
            int randWord = Gibbs.getRandWordInd(lWordsPr);
            // replacing the word in S (10.9 7)
            S.set(randSeq, (String)lWords.get(randWord));
            
            // Generate new motif and background (10.9 7)
            Gibbs.getTheta();
            Gibbs.getThetaZero();
            
            // Update overall score
            score2 = Gibbs.logScoreSumS();
            if((int)score == (int)score2) {
                counter++;
            } else {
                counter = 0;
            }
            score = score2;
            
            // End the loop when score does not change Figure 10.9 8
            if(counter==20){break;} // if the score the same a twenty times, stop
        }
        
        /* Optional debug printouts */
        System.out.println("--------------------------------------------------");
        Gibbs.printTheta();
        System.out.println("--------------------------------------------------");
        Gibbs.printThetaZero();
        System.out.println("--------------------------------------------------");
        Gibbs.printSmotif();
        System.out.println("--------------------------------------------------");
        Gibbs.printS();
        System.out.println("--------------------------------------------------");
        Gibbs.printMotifOnly();
        System.out.println("--------------------------------------------------");
        Gibbs.printConcensus();
        System.out.println("--------------------------------------------------");
        
        // create output path and write to file
        out = Gibbs.promptOutFile();
        try (PrintWriter outputStream = Gibbs.outputStreamer(out.getName())) {
            Gibbs.writeMotifToFile(outputStream);
        }
    }  
}