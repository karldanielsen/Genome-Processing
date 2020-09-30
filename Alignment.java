/*
 *  This program recieves a reference genome as a file argument,
 *  and then accepts reads from another, unsequenced genome (from 
 *  STDIN)and locates places where the two differ.
 *  
 *  Multiple differences are checked for:
 *  SNP (Snips)       where one allele is changed
 *  INS (Insertions)  where a sequence is inserted in the new genome
 *  DEL (Deletions)   where a sequence is missing in the new genome
*/

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class Alignment{

    private final int K = 50;
    private final int MERLEN = 17;
    private final int MAPSIZE = 1000000;
    
    private HashMap<String,Integer> ref;
    private HashMap<String,Integer> kmers;
    private HashMap<String,Integer> hasMatch;

    private ArrayList<String> unmatcheds;

    private String genome;
    private String genomeName;
    
    public static void main(String[] args){
	ref   = new HashMap<String,Integer>(MAPSIZE/3);
	kmers = new HashMap<String,Integer>(MAPSIZE);
	hasMatch = new HashMap<String,Integer>(MAPSIZE);
	genome = "";

	//Read in the reference genome
	Scanner scanner = null;
	try{
	    scanner = new Scanner(
			      new File(args[0]));
	}
	catch(FileNotFoundException e){
	    System.out.println("File not found.");
	}

	//Prep reference genome for easy comparison to new genome
	populateRefs(scanner);

	//Attempt to match reads to locations in the reference genome
	unmatcheds = indexReads();

	//Now, check unmatched reads for each type of error
	ArrayList<String> snpOutput = new ArrayList<String>();
	ArrayList<String> delOutput = new ArrayList<String>();
	ArrayList<String> insOutput = new ArrayList<String>();
	ArrayList<String> testOutput = new ArrayList<String>();
	
	String[] submers;
	boolean foundSNP = false;
	for(String unmatched : unmatcheds){
	    //Divide the unmatched read into 3 mers of length 17.
	    try{
		submers = new String[]{unmatched.substring(0,MERLEN),
				       unmatched.substring(MERLEN,MERLEN*2),
				       unmatched.substring(MERLEN*2-1,K)};
	    }
	    catch(StringIndexOutOfBoundsException e){
	        continue;
	    }

	    //index vars track where in the reference genome the mers match
	    int indexOfMer1 = -1;
	    int indexOfMer2 = -1;
	    int indexOfMer3 = -1;
	    int mersFound = 0;

	    for(int i = 0; i < submers.length; i++){
		//If two mers of a read match but one doesn't,
		//check it for problems.
		if(kmers.containsKey(submers[i])){
		    mersFound++;
		    if(i == 0)
			indexOfMer1 = kmers.get(submers[i]);
		    else if(i == 1)
			indexOfMer2 = kmers.get(submers[i])-MERLEN;
		    else
		    indexOfMer3 = kmers.get(submers[i])-(MERLEN*2-1);
		}
	    }
	    
	    //If two of the kmers match, check for SNPs
	    if(mersFound > 1)
		foundSNP = checkSNP(snpOutput,unmatched,indexOfMer1,indexOfMer2,indexOfMer3);

	    //If no SNP is found, check for INS/DEL
	    if(!foundSNP)
		checkIndel(insOutput,delOutput,unmatched,indexOfMer1);
	}

	/*
	 *  Prints according to FASTA formatting instructions here
	 *     --> https://cm122.herokuapp.com/ans_file_doc
	*/
	System.out.println(genomeName);
	
	//Only print SNPs with 3+ appearances
	System.out.println(">SNP");
	auditErrors(snpOutput,2);
	
	//Only print insertions with 2+ appearances
	System.out.println(">INS");
	auditErrors(insOutput,1);

	//Only print deletions with 2+ appearances
	System.out.println(">DEL");
	auditErrors(delOutput,1);
	
    }

    /* 
     *  Steps through the file, and populates
     *  "global" variables with the full sequence, read-length,
     *  and k-mer length substrings of the reference genome.
     */
    private void populateRefs(Scanner scanner){
	int refIndex  = 0;
	int kmerIndex = 0;
	String line;
	String previous = "";
	scanner.nextLine();
	while(scanner.hasNextLine()){
	    line = scanner.nextLine();
	    genome += line;

	    //Populate reference genome reads
	    if(!previous.equals(""))
		for(int i = line.length()-K; i < line.length(); i++){
		    ref.putIfAbsent(previous.substring(i,line.length()) +
				    line.substring(0,i+K-line.length()),
				    refIndex++);
		}
	    for(int i = 0; i < line.length()-K; i++){
		ref.putIfAbsent(line.substring(i,i+K),refIndex++);
	    }

	    //Populate reference genome k-mers (length 17)
	    //Also in charge of incrementing refIndex
	    if(!previous.equals(""))
		for(int i = line.length()-MERLEN; i < line.length(); i++){
		    kmers.putIfAbsent(previous.substring(i,line.length()) +
				      line.substring(0,i+MERLEN-line.length()),
				      kmerIndex++);
		}
	    for(int i = 0; i < line.length()-MERLEN; i++)
		kmers.putIfAbsent(line.substring(i,i+MERLEN),kmerIndex++);
	    previous = line;
	}
	scanner.close();
    }


    /*
     *  Compares each new read to the hashmap of 
     *  reference reads. If the read is in the hashmap, it can 
     *  be ignored. If not, it is returned so the difference can be analysed.
     */
    private ArrayList<String> indexReads(){
	Scanner scanner = new Scanner(System.in);
	String read = "";
	String pair = "";
	ArrayList<string > unmatched = new ArrayList<String>();

	genomeName = scanner.nextLine();
	
	//For each read,attempt to index it into the hashmap
	while(scanner.hasNextLine()){
	    String[] lines = scanner.nextLine().split(",");
	    read = lines[0];
	    if(lines.length < 2){
		if(!ref.containsKey(read))
		    unmatched.add(read);
		continue;
	    }
	    pair = lines[1];
	    //Unmatched reads must be checked for errors.
	    //If they match an error time, a correctly formatted message
	    //is printed, otherwise, it is treated like garbage and thrown out
	    if(!ref.containsKey(read))
		unmatched.add(read);
	    if(!ref.containsKey(pair))
	    	unmatched.add(pair);
	}
	scanner.close();
	return unmatched;
    }

    /*
     *  Checks a given unmatched read for places substrings of it (mers)
     *  match the reference genome. If 2/3 mers match, the differences are
     *  added to the list of SNPs, and true is returned. Otherwise, false
     *  is returned.
     */
    private boolean checkSNP(ArrayList<String> snpOutput. String unmatched, int indexOfMer1, int indexOfMer2, int indexOfMer3){
	//First check the three k-mers of the read
	int indexOfRead = -1;
	if((indexOfMer1 == indexOfMer2 ||
	    indexOfMer2 == indexOfMer3 ||
	    indexOfMer3 == indexOfMer1)){
	    if(indexOfMer1 > -1)
		indexOfRead = indexOfMer1;
	    else
		indexOfRead = indexOfMer2;
	}
	//If a match to the unmatched read is found, check for the error
	if(indexOfRead >= 0){
	    String OGRead = genome.substring(indexOfRead,indexOfRead+K);
	    
	    //First, check for SNP's
	    String snpout = "";
	    int snpfound = 0;
	    for(int i = 0; i < K; i++){
		if(OGRead.charAt(i) != unmatched.charAt(i)){
		    //Must check the case where multiple non-consecutive errors are found
		    if(snpfound > 2){

			testOutput.clear();
			break;
		    }
		    testOutput.add(OGRead.charAt(i) + "," + unmatched.charAt(i) + "," + (indexOfRead+i));
		    snpfound++;
		}
	    }
	    for(String to : testOutput)
		snpOutput.add(to);
	    testOutput.clear();
	    return true;
	}
	return false;
    }
    /*
     *  Checks a given unmatched read for places where there is a gap
     *  in either the reference or new genome, eg.
     *  ABCDEFGHIJK
     *  ABCD___EFGH
     *
     *  If a gap is found, it is added to the INS or DEL data structures.
     */
    private void checkIndel(ArrayList<String> insOutput, ArrayList<String> delOutput, String unmatched, int indexOfMer1){
	boolean legalChange = true;
	for(int i = 0; i < K; i++){
	    String OGRead = genome.substring(indexOfMer1,indexOfMer1+K);
	    if(OGRead.charAt(i) != unmatched.charAt(i)){
		//If a mismatch is found, check every index after it, to see if the sequence repeats

		//This finds insertions, where unmatched+c = OGRead for all indices past the mismatch
		for(int c = 0; c < 5; c++){
		    legalChange = true;
		    if(i+c > K)
			break;
		    for(int j = i; j < K-c; j++){
			if(OGRead.charAt(j) != unmatched.charAt(j+c))
			    legalChange = false;
		    }
		    if(legalChange)
			insOutput.add(unmatched.substring(i,i+c) + "," + (indexOfMer1+i)); 
		}
			
		//This finds deletions, where unmatched = OGRead+c for all indices past the mismatch
		for(int c = 0; c < 5; c++){
		    legalChange = true;
		    if(i+c > K)
			break;
		    for(int j = i; j < K-c; j++){
			if(OGRead.charAt(j+c) != unmatched.charAt(c))
			    legalChange = false;
		    }
		    if(legalChange)
			delOutput.add(OGRead.substring(i,i+c) + "," + (indexOfMer1+i)); 
		}
		break;
	    }
	}
    }
    
    /*
     *  Some read errors may be due to tool malfunction. This method
     *  only prints errors with > min appearances in the check methods.
     */
    private void auditErrors(ArrayList<String> output, int min){
	HashMap<String,Integer> appearances = new HashMap<String,Integer>(10000);
	HashMap<String,Boolean> printed = new HashMap<String,Boolean>(10000);

	//Count an error's detected appearances
	for(String out : output){
	    if(appearances.get(out) != null)
		appearances.put(out,appearances.get(out)+1);
	    else
		appearances.put(out,0);
	}
	for(String out : output){
	    if(printed.get(out) != null)
		continue;
	    if(appearances.get(out) > min){
		System.out.println(out);
		printed.put(out,true);
	    }
	}
    }

    /*
     *  Clear the data structures that contain information on the ref genome.
     */
    private void clearRefGenome(){
	ref.clear();
	kmers.clear();
	hasMatch.clear();
	genome = "";
    }

    /*
     *  Clear the data structures that contain information on the new genome.
     */
    private void clearNewGenome(){
	unmatcheds.clear();
    }
}
