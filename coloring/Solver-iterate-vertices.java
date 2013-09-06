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
    //Stack<Integer> vStack = new Stack<Integer>();
    Queue<Integer> vQueue = new Queue<Integer>();
    
    vArray = new Vertix[g.V()];
    for (int v = 0; v < g.V(); v++) {
      vArray[v] = new Vertix(v);
    }
    java.util.Arrays.sort(vArray);
    
    for (int k = 0; k < g.V(); k++) {
      // find the vertix to color
      if (vArray[k].color() != -1) continue;
      
      //vStack.push(k);
      vQueue.enqueue(k);
      
      //while (vStack.size() > 0) {
      while (vQueue.size() > 0) {
        //int v = vStack.pop();
        int v = vQueue.dequeue();
            
        // try to color v with an existing color
        int color = -1;
        for (int c = 0; c < numUsedColors; c++) {
          if (colorVertix(v, c)) {
            color = c;
            break;
          }
        }
      
        // if an existing color could not be used, increment a new color and use it
        if (color == -1) {
          numUsedColors++;
          color = numUsedColors - 1;
          colorVertix(v, color);      
        }
      
      
        // color as many uncolored vertices as possible
        for (int i = 0; i < g.V(); i++) {
          if (i == v) continue;
          if (vArray[i].color() == -1) {
            colorVertix(i, color);
            //vStack.push(i);
            vQueue.enqueue(i);
          }
        }
      }
    }
        
    //printVColor();

    // prepare the solution in the specified output format        
    System.out.println(numUsedColors + " 0");
    String s = "";
    int[] outColor = new int[g.V()];
    for (int v = 0; v < g.V(); v++) {
      outColor[vArray[v].vertix()] = vArray[v].color();
    }
    for (int v = 0; v < g.V(); v++) {
      s+= outColor[v] + " ";
    }
    System.out.println(s);
    
    printViolations();

  }
     
  private boolean colorVertix(int v, int c) {
    // colors vertix v with color c if possible (violates no constraint)
    // returns true for success, false for failure
    // note: will recolor if v is already colored

    for (int w : g.adj(v)) if (vArray[w].color() == c) return false;  // not feasible
    vArray[v].setColor(c);
    numVColored++;
    if (!usedColors[c]) usedColors[c] = true;
    return true;     
  }
  
  private void printViolations() {
    boolean flag = false;
    String s = "";
    for (int v = 0; v < g.V(); v++) {
      for (int w : g.adj(v)) {
        if (vArray[w].color() == vArray[v].color()) {
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
    int color;
    public Vertix(int v) {
      vertix = v;
      color = -1;
      dimension = 0;
      for (int i : g.adj(v)) {
        dimension++;
      }
    }
    public int compareTo(Vertix that) {
      // reverse order
      if (this.dimension() > that.dimension()) return -1;
      if (this.dimension() < that.dimension()) return 1;
      return 0;
    }
    public int vertix() { return vertix; }
    public int dimension() { return dimension; }
    public void setColor(int c) { color = c; }
    public int color() { return color; }
    public String toString() {
      return "vertix, dimension, color: " + vertix + ", " + dimension + ", " + color;
    }
  }
  

}