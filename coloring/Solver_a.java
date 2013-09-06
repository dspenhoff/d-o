import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.BitSet;

/**
 * The class <code>Solver</code> is an implementation of a greedy algorithm to solve the knapsack problem.
 *
 */
public class Solver_a {
  
  int N, E;
  Graph g;
  int numUsedColors;
  int numVColored;
  int[] vColor;
    
    /**
     * The main class
     */
    public static void main(String[] args) {
        try {
          Solver solver = new Solver(args);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Read the instance, solve it, and print the solution in the standard output
     */
    public Solver_a(String[] args) throws IOException {
        String fileName = null;
        
        // get the temp file name
        for(String arg : args){
            if(arg.startsWith("-file=")){
                fileName = arg.substring(6);
            } 
        }
        if(fileName == null)
            return;
        
        // read the lines out of the file
        List<String> lines = new ArrayList<String>();

        BufferedReader input =  new BufferedReader(new FileReader(fileName));
        try {
            String line = null;
            while (( line = input.readLine()) != null){
                lines.add(line);
            }
        }
        finally {
            input.close();
        }
              
        // parse the data in the file
        String[] firstLine = lines.get(0).split("\\s+");
        N = Integer.parseInt(firstLine[0]);
        E = Integer.parseInt(firstLine[1]);
        g = new Graph(N);

        for(int i= 0; i < E; i++){
          String line = lines.get(i+1);
          String[] parts = line.split("\\s+");
          int v = Integer.parseInt(parts[0]);
          int w = Integer.parseInt(parts[1]);
          g.addEdge(v, w);
        }
        
        //System.out.println(g);
        
        // note: g.V() returns # vertices; g.E() returns # edges; vertix.adj() returns an iterable object
        // containing vertices adjacent to vertix
        
        // color the vertices
        // todo
        /*
        idea for a greedy heuristic:
        - call the colors 0, 1, 2, ...
        - array vertixColors has the color for each vertix that has been colored, 0 if uncolored
        - start with vertix 0 and assign it color 0 . push v0 on a stack. (alternately, try a queue or random selection)
        - while the stack is not empty, pop a vertix (alternately dequeue, e.g., test LIFO and FIFO)
        - get v.adj(); 
        */
        
        vColor = new int[g.V()];
        for (int v = 0; v < g.V(); v++) vColor[v] = -1;
        numVColored = 0; 
        numUsedColors = 0;
        boolean[] usedColors = new boolean[1000];  
      
        while (numVColored < g.V()) { // control loop for disconnected graphs
        //System.out.println("top subgraph control loop, numVColored = " + numVColored);
          
          int startingV = 0;
          for (int i = 0; i < g.V(); i++) {
            if (vColor[i] == -1) {
              startingV = i;
              break;
            }
          }
          vColor[startingV] = 0;
          numVColored++;
          if (!usedColors[0]) {
            usedColors[0] = true;
            numUsedColors++;
          }
          //Stack<Integer> vNext = new Stack<Integer>();        
          //vStack.push(startingV);
          Queue<Integer> vNext = new Queue<Integer>();
          vNext.enqueue(startingV);
          
          while (vNext.size() > 0) {  // control loop for a connected subgraph
                    
            //int v = vNext.pop();
            int v = vNext.dequeue();
                   
            // what colors have been used by adjacent nodes?
            usedColors[vColor[v]] = true;
            for (int i : g.adj(v)) { 
              if (vColor[i] >= 0) {
                if (vColor[i] == vColor[v]) {
                  vColor[i] = -1;
                  numVColored--;
                  break;
                }
 
                usedColors[vColor[i]] = true; 
              }
            }
          
            // determine the color to use
            int theColor = -1;
            for (int c = 0; c < numUsedColors; c++) {
              if (!usedColors[c]) {
                theColor = c;
                break;
              }
            }
            if (theColor == -1) {
              numUsedColors++;
              theColor = numUsedColors - 1;
            }
            
            // apply the color to uncolored adjacent nodes
            for (int w : g.adj(v)) {
              if (vColor[w] == -1) {
                vColor[w] = theColor;
                //System.out.println("colored vertix = " + w + " from vertix = " + v);
                numVColored++;
                //vNext.push(w);
                vNext.enqueue(w);
              }
            }
            //printVColor(); 
            
          }
          //System.out.println("bottom subgraph control loop, numVColored = " + numVColored);      
        }
        
        printVColor();
  
      
        // prepare the solution in the specified output format
        // todo
        
        System.out.println(numUsedColors + " 0");
        String s = "";
        for (int v = 0; v < g.V(); v++) {
          s += vColor[v] + " ";
        }
        System.out.println(s);
        
        printViolations();

    }
    
    private void printVColor() {
      System.out.println("num colors = " + numUsedColors);
      for (int v = 0; v < g.V(); v++) {
        System.out.println("vertix, color = " + v + ", " + vColor[v]);
      }
    }
    
    private void printViolations() {
      String s = "";
      for (int v = 0; v < g.V(); v++) {
        for (int w : g.adj(v)) {
          if (vColor[w] == vColor[v]) {
            s += "(" + v + "," + w + ")";
          }
        }
      }
      System.out.println("Constraint violations: " + s);
    }

}