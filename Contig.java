/*
 *  This program tries to construct a full genome from reads.
 *  It recieves FASTA-style input and converts it into a graph input 
 *  format, which looks like:
 *  1 -> 2,3,4
 *  2 -> 1,4
 *  with numbers representing nodes (reads) and positioning representing
 *  directed edges. A read points to another read if they overlap, eg.
 *  ABCDE would point to BCDEF, as their suffixes and prefixes overlap.
 *
 *  After the input has been formatted, a Debrujin analysis is done  by
 *  finding every path that only goes through 1-1 nodes from non 1-1 nodes.
 *  Finally, isolated cycles are found by checking 1-1 nodes which have not 
 *  been visited by prior analysis.
 *
 *  The resulting FASTA-formatted data file contains every confidently
 *  reconstructed section of the genome.
 */

import java.util.*;
public class Contig{
    private static ArrayList<String[]> graph;
    private static boolean[] visited;
    private static final int NODES = 20000;
    private static final int THRESHOLD = 2;
    private static final int MERLEN = 30;
    private static String[] mers;
	
    public static void main(String[] args){
	reads = new String[NODES];
	graph = new ArrayList<String[]>();
	visited = new boolean[NODES];

	String graphAsString = "";
	int[] edgesTo = new int[NODES];
	int[] edgesFrom = new int[NODES];
	Scanner scanner = new Scanner(System.in);
	
	//Code to sort input into a graph
	String firstLine = scanner.nextLine().split(",")[0];
	int K = firstLine.length();
	reads[0] = firstLine;
	for(int i = 1; scanner.hasNextLine(); i++){
	    reads[i] = scanner.nextLine().split(",")[0];
	}

	//Preprocess reads by removing those with low-appearance k-mers.
	//This can be done by populating a hashtable of kmers-># appearances,
	//And then looping through it and deleting reads with low appearances
	Hashtable<String,Integer> kmers = new Hashtable<String,Integer>(NODES);
	for(String read : reads){
	    if(read.length() < 50)
		    continue;
	    for(int i = 0; i <= (K - MERLEN); i++){
		if(kmers.get(read.substring(i,i+MERLEN)) == null)
		    kmers.put(read.substring(i,i+MERLEN),1);
		else
		    kmers.put(read.substring(i,i+MERLEN),kmers.get(mer.substring(i,i+MERLEN))+1);
					
	    }
	}

	//With the hashtable populated, it's time to remove reads
	//I remove them by setting their value to an illegal string
	for(int i = 0; i < reads.length; i++){ //cannot be done w/ for:each b/c item is mutated
	    if(reads[i].length() < 50){
		reads[i] = "OFF";
		continue;
	    }
	    for(int j = 0; j <= (K-MERLEN); j++)
	    	if(kmers.get(reads[i].substring(j,j+MERLEN)) < THRESHOLD){
	    	    reads[i] = "OFF"; break;}
	}

	//Remove duplicate reads. A genome sample will obviously have
	//multiple copies of the same section of RNA. Only keep one.
	for(int suf = 0; suf < reads.length; suf++){
	    if(reads[suf].equals("OFF"))
		continue;
	    for(int pre = suf+1; pre < reads.length; pre++){
		if(reads[suf].equals(reads[pre])){
			reads[pre] = "OFF";
		    }
	    }
	}
		
	//Write the graph as a giant string in the graph-format specified.
	for(int suf = 0; suf < reads.length; suf++){
	    if(reads[suf].equals("OFF"))
		continue;
	    boolean matchFound = false;
	    graphAsString += suf + " ->";
	    for(int pre = 0; pre < reads.length; pre++){
		if(reads[pre].equals("OFF"))
		    continue;
		if(reads[suf].substring(1,K).equals(reads[pre].substring(0,K-1))){
		    graphAsString += " " + pre;
		}
	    }
	    
	    graphAsString += "\n";
	}
	scanner.close();

	
	//Populate edgesTo(prefix)[] edgesFrom(prefix)[]
	//It is important to know how many reads point to and from
	//nodes to find 1-1 nodes and run DeBrujin
	for(int j = 0; j < reads.length; j++){
	    if(reads[j].equals("OFF"))
		continue;
	    for(int k = 0; k < reads.length; k++){
		if(reads[k].equals("OFF"))
		    continue;
		if(reads[k].substring(0,K-1).equals(reads[j].substring(0,K-1)))
		    edgesFrom[k] += 1;
		if(reads[k].substring(0,K-1).equals(reads[j].substring(1,K)))
		    edgesTo[k] += 1;
	    }
	}

	
	scanner = new Scanner(graphAsString);
	
	//Code to parse out graph from graph-formatted input
	while(scanner.hasNextLine()){
	    String[] next = scanner.nextLine().split("(,| )");
	    graph.add(next);
	}
	Collections.sort(graph, new SortFirstIndex());
       	scanner.close();


	//DeBrujin code
	//Analysis done according to online pseudocode,
	//Traverses forward from non 1-1 nodes until another
	//non 1-1 node is found
	ArrayList<String> paths = new ArrayList<String>();
	for(String[] node : graph){
	    if(edgesFrom[Integer.parseInt(node[0])] != 1 ||
	       edgesTo[Integer.parseInt(node[0])] != 1){
		visited[Integer.parseInt(node[0])] = true;
		if(node.length == 2){		    
		    paths.add(node[0]);
		    continue;
		}    
		    if(edgesFrom[Integer.parseInt(node[0])] > 0){
		    for(int i = 2; i < node.length; i++){
			String newNode = node[i];
			String newPath = node[0];
			while(edgesFrom[Integer.parseInt(newNode)] == 1 &&
			      edgesTo[Integer.parseInt(newNode)] == 1){
			    visited[Integer.parseInt(newNode)] = true;
			    newPath+= " -> " + newNode;
			    if(getNode(newNode).length != 2)
				newNode = getNode(newNode)[2];
			    else{ //This secondary case is important, it catches single-node paths
				newNode = reads[Integer.parseInt(newNode)].substring(K-1,K);
				break;
			    }
			
			}
			try{
			    visited[Integer.parseInt(newNode)] = true;
			}
			catch(NumberFormatException e){}
			newPath += " -> " + newNode;
			paths.add(newPath);
		    }
		    
		}
	    }
	}	
	ArrayList<String> output = new ArrayList<String>();
	ArrayList<String> used = new ArrayList<String>();
	int previous = 0; //assume starting at reads[0]

	//Some isolated cycles are not findable by traversing from the
	//graph's start, these are added separately
	paths = addIsolated(paths);

	//All this simply formats the output back to FASTA
	String newOut = "";
	for(String path : paths){
	    String[] out = path.split(" ");
	    if(!(previous == Integer.parseInt(out[0])))
		used = new ArrayList<String>();
	    newOut+=(reads[Integer.parseInt(out[0])].substring(0,K-1));
	    int i = 2;
	    for(; i < out.length;i+=2){
		try{
		newOut+=(reads[Integer.parseInt(out[i])].substring(K-2,K-1));
		}
		catch(NumberFormatException e){
		    newOut+=out[i];
		}
		}
	    if(out.length == 1)
		newOut = (reads[Integer.parseInt(out[0])]);
	    if(!used.contains(newOut))
		output.add(newOut);
	    previous = Integer.parseInt(out[0]);
	    used.add(newOut);
	    newOut="";
	}

	//Add in isolated nodes that did not have the original
	//path detection run for them, as they
	//have no edges leaving them
	for(int i = 0; i < edgesFrom.length;i++){
	    if(edgesFrom[i] == 0 && !reads[i].equals("OFF"))
		output.add(reads[i]);
	}
	//Finally, print all nodes. The chromosome name and word
	//"ASSEMBLY" are added in post-processing
	for(String out : output)
	  System.out.println(out);
    }


    /*
     *  In order to find self-contained cycles, where every node has
     *  one edge entering it and leaving it, nodes not visited by the
     *  DeBrujin analysis are checked for cycles.
     */
    public static ArrayList<String> addIsolated(ArrayList<String> paths){
	//Find an unvisited node in graphs
	for(String[] node : graph){
	    if(!visited[Integer.parseInt(node[0])]){
		visited[Integer.parseInt(node[0])] = true;
		String newNode = node[2];
		String newPath = node[0] + " -> " + newNode;
		//Travel forward until an already encountered node is found
		while(!visited[Integer.parseInt(newNode)]){
		    visited[Integer.parseInt(newNode)] = true;
		    newPath+= " -> " + getNode(newNode)[2];
		    newNode = getNode(newNode)[2];
		}
		paths.add(newPath);
	    }
	}
	return paths;
    }
    public static String[] getNode(String num){
	for(String[] node : graph){
	    if(num.equals(node[0]))
		return node;
	}
	return new String[]{"1","2","3","4","5"};
    }
}
