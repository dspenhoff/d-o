import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.BitSet;

/**
 * The class <code>Solver</code> is an implementation of a greedy algorithm to solve the knapsack problem.
 *
*/

public class Solver {
  
  int N, E;
  Graph g;
  int numUsedColors;
  int numVColored;
  int[] vColor;
  int MAX_COLORS = 1000;
  boolean[] usedColors = new boolean[MAX_COLORS];  
  Vertix[] vArray;
    
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
  public Solver(String[] args) throws IOException {
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
        
    vColor = new int[g.V()];
    for (int v = 0; v < g.V(); v++) vColor[v] = -1;
    numUsedColors = 0;
  
    /*
      this version simply iterates through the vertices
      other possible versions:
      - sort the vertices by number of adjacent vertices
      - push each colored vertix on a stack and pop through the vertices 
      - use a queue instead of a stack
      - use a randomizing queue
      - (need to account for disconnected subgraphs)
      regardless, the 1-2-3 steps below for coloring a given vertix should be the same
    */
    

    vArray = new Vertix[g.V()];
    for (int v = 0; v < g.V(); v++) {
      vArray[v] = new Vertix(v);
    }
    java.util.Arrays.sort(vArray);
    
    for (int i = 0; i < g.V(); i++) {
      colorVertix(vArray[i].vertix());      
    }
    
    /*
    for (int v = 0; v < g.V(); v++) {
      colorVertix(v);      
    }
    */
        
    //printVColor();

    // prepare the solution in the specified output format        
    System.out.println(numUsedColors + " 0");
    String s = "";
    for (int v = 0; v < g.V(); v++) {
      s += vColor[v] + " ";
    }
    System.out.println(s);
    
    printViolations();

  }
    
  private void colorVertix(int v) {
    // color vertix v:
    // 1. determine the minumum unused color for vertices adjacent to v: c (this is most of the work)
    // 2. color v with c
    // 3. update counters, etc.
    
    int colorForV = -1;    // default color
    
    boolean[] adjUsedColors = new boolean[numUsedColors];
    for (int w : g.adj(v)) {
      if (vColor[w] != -1) {
        adjUsedColors[vColor[w]] = true;
      }
    }
    for (int c = 0; c < numUsedColors; c++) {
      if (!adjUsedColors[c]) {
        colorForV = c;
        break;
      }
    }
    if (colorForV < 0) {
      colorForV = numUsedColors++;
    }
    
    vColor[v] = colorForV;
    if (!usedColors[colorForV]) {
      usedColors[colorForV] = true;
    }      
  }

  private void printVColor() {
    System.out.println("num colors = " + numUsedColors);
    for (int v = 0; v < g.V(); v++) {
      System.out.println("vertix, color = " + v + ", " + vColor[v]);
    }
  }
    
  private void printViolations() {
    boolean flag = false;
    String s = "";
    for (int v = 0; v < g.V(); v++) {
      for (int w : g.adj(v)) {
        if (vColor[w] == vColor[v]) {
          flag = true;
          s += "(" + v + "," + w + ")";
        }
      }
    }
    if (flag) System.out.println("Constraint violations: " + s);
  }
  
  private class Vertix implements Comparable<Vertix> {
    private int vertix;
    private int dimension;
    public Vertix(int v) {
      vertix = v;
      dimension = 0;
      for (int i : g.adj(v)) {
        dimension++;
      }
    }
    public int compareTo(Vertix that) {
      // reverse order
      if (this.vertix() > that.vertix()) return -1;
      if (this.vertix() < that.vertix()) return 1;
      return 0;
    }
    public int vertix() { return vertix; }
    public int dimension() { return dimension; }
  }
  

}