/**
 * Basic implementation of the DeBrujin algorithm, later modified for
 * use in Contigs.java. 
 * <p>
 * Finds paths by investigating every node that isn't a through node
 * (one input one output). Then finds isolated cycles by investigating
 * unvisited nodes.
 *
 * @author Karl Danielsen
 * @version 0.1
 */

import java.util.*;

class DeBrujin{
    private static ArrayList<String[]> graph;
    private static final int NODES = 512;

    /**
     * Performs the analysis using stdin as input.
     *
     * @param args passed-in command line arguments.
     */
    public static void main(String[] args){
        String[] mers = new String[NODES];
        String graphAsString = "";
        graph = new ArrayList<String[]>();
        Scanner scanner = new Scanner(System.in);
        int[] edgesFrom = new int[NODES];
        int[] edgesTo = new int[NODES];

        //PreProcessing to get debrujin data in euler path form
        int K = Integer.parseInt(scanner.nextLine());
        for(int i = 0; scanner.hasNextLine(); i++){
            mers[i] = scanner.nextLine();
        }
        for(int suf = 0; suf < mers.length; suf++){
            boolean matchFound = false;
            graphAsString += suf + " ->";
            for(int pre = 0; pre < mers.length; pre++){
                if(mers[suf].substring(1,K).equals(mers[pre].substring(0,K-1))){
                    matchFound = true;
                    graphAsString += " " + pre;
                }
            }
            if(matchFound == false)
                graphAsString += " NIL";
            graphAsString += "\n";
        }
        scanner.close();
        scanner = new Scanner(graphAsString);
        //Original Euler Path Code
        while(scanner.hasNextLine()){
            String[] next = scanner.nextLine().split("(,| )");
            if(next[2].equals("NIL"))
                edgesFrom[Integer.parseInt(next[0])] = 0;
            else
                edgesFrom[Integer.parseInt(next[0])] = next.length-2;
            for(int i = 2; i < next.length; i++){
                if(!next[i].equals("NIL"))
                    edgesTo[Integer.parseInt(next[i])]+=1;
            }
            graph.add(next);
        }

        int start=0;
        int end=0;
        for(int i = 0; i < edgesFrom.length; i++){
            if(edgesFrom[i] > edgesTo[i])
                start = i;
            if(edgesFrom[i] < edgesTo[i])
                end = i;
        }
        Collections.sort(graph, new SortFirstIndex());
        for(String [] eh : graph){
            for(String meh : eh)
                System.out.print(meh + " ");
            System.out.println();
        }
        scanner.close();

        //Now have an arraylist where each index corresponds to an array of edges

        //Find the first cycle
        System.out.println(start + " " + end);
        ArrayList<ArrayList<String>> visited = new ArrayList<ArrayList<String>>();
        for(int i = 0; i < graph.size(); i++)
            visited.add(new ArrayList<String>());
        String output = graph.get(start)[0] + " -> "+ findCycle(Integer.parseInt(graph.get(start)[2]),
                                                                end, visited);
        graph.get(start)[2] = "NIL";
        ArrayList<String> cycle = new ArrayList<String>();
        String[] cyc = output.split(" -> ");
        Collections.addAll(cycle,cyc);
        //System.out.println(findCycle(8,6));

        //Repeat until every edge has been used
        while(nonNilExists()){
            for(int i = 0; i < cycle.size(); i++){
                int spot = Integer.parseInt(cycle.get(i));
                for(int j = 2; j < graph.get(spot).length; j++){
                    if(graph.get(spot)[j].equals("NIL"))
                        continue;
                    visited = new ArrayList<ArrayList<String>>();
                    for(int f = 0; f < graph.size(); f++)
                        visited.add(new ArrayList<String>());
                    String newCycle = spot + " -> " + findCycle(Integer.parseInt(graph.get(spot)[j]),spot,visited);
                    graph.get(spot)[j] = "NIL";
                    if(newCycle.equals(""))
                        continue;
                    String[] newCyc = newCycle.split(" -> ");
                    for(int k = newCyc.length-2; k > -1; k-- )
                        cycle.add(i,newCyc[k]);
                }
            }
        }
        System.out.print(mers[Integer.parseInt(cycle.get(0))]);
        for(int i = 1; i < cycle.size();i++)
            System.out.print(mers[Integer.parseInt(cycle.get(i))].charAt(K-1));
        System.out.println();
    }

    /** 
     * Finds a cycle from the given start index. Called for every start index,
     * and then for every unvisited node.
     *
     * @params start the node that is currently being investigated. Updates 
     *               with recursive calls.
     * @params goal the node to return to. Equal to initial call's "start".
     * @params visited the list of visited nodes, to be avoided in searches.
     */
    public static String findCycle(int start, int goal, ArrayList<ArrayList<String>> visited){
        if(goal == start)
            return " -> " + goal;
        for(int i = 2; i < graph.get(start).length; i++){ //Check each edge
            if(graph.get(start)[i].equals("NIL"))
                continue;
            if(visited.get(start).contains(graph.get(start)[i]))
                continue;
            visited.get(start).add(graph.get(start)[i]);
            String currCycle = findCycle(Integer.parseInt(graph.get(start)[i]),goal,visited);
            if(Integer.parseInt(graph.get(start)[i]) == goal){
                String temp = graph.get(start)[i];
                graph.get(start)[i] = "NIL";
                return start + " -> " + temp;
            }
            else if(currCycle.contains(Integer.toString(goal))){
                String temp = graph.get(start)[i];
                graph.get(start)[i] = "NIL";
                return start + " -> " + currCycle;
            }
        }
        return "";
    }
    /**
     * Checks whether their are still nodes not in cycles within the graph.
     * 
     * @return true if the program must continue searching, false otherwise.
     */
    public static boolean nonNilExists(){
        for(int i = 0; i < graph.size(); i++)
            for(int j = 2; j < graph.get(i).length; j++)
                if(!graph.get(i)[j].equals("NIL"))
                    return true;
        return false;
    }
}

/**
 * Comparator for arrays of strings, used to topographically sort data.
 */
class SortFirstIndex implements Comparator<String[]>{
    public int compare(String[] a, String[] b){
        if(Integer.parseInt(a[0]) > Integer.parseInt(b[0]))
            return 1;
        return -1;
    }
}
